import os
import pandas as pd
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import make_scorer, recall_score, confusion_matrix, classification_report
from joblib import parallel_backend
# to do further optimise features to choose.

class PeakListClassifier():
    def __init__(self, input_file):
        directory = os.path.dirname(input_file)
        df = pd.read_csv(input_file)

        # define features and target
        features = [
            "charge",
            "fit_score",
            "interference",
            "log_s2n",
            "pc_missing_peaks",
            "mass_error_std",
            "chisq_p_value",
            "correlation_p_value",
            "correlation_coefficient",
            "chisq_stat"
        ]
        target = "validated"

        df = df.drop(df.index[-1]) # remove last line - no entry
        checked_data = df[df["checked"] == True]
        other_data = df[df["checked"] != True]
        combined_data = pd.concat([checked_data, other_data], ignore_index=True)

        # extract feature and target values for all data, then apply scaling
        x_combined = combined_data[features]
        y_combined = combined_data[target]
        scaler = StandardScaler()
        x_combined_scaled = scaler.fit_transform(x_combined)

        # divide data up again after scaling
        x_checked_scaled = x_combined_scaled[:len(checked_data)]
        y = y_combined[:len(checked_data)] # training/test set data

        # split data 75/25
        x_train_scaled, x_test_scaled, y_train, y_test = train_test_split(
            x_checked_scaled,
            y,
            test_size=0.25
        )

        y_train = y_train.astype('bool')
        y_test = y_test.astype('bool')

        # train Logistic Regression Model
        logistic_model = self.train_logistic_regression(x_train_scaled, y_train)
        joblib.dump(logistic_model, os.path.join(directory, "logisticRegressionModel.pk1"))

        # train Gradient Boosting Model
        gradient_boosting_model = self.train_gradient_boosting(x_train_scaled, y_train)
        joblib.dump(gradient_boosting_model, os.path.join(directory, "gradientBoostingModel.pk1"))

        self.evaluate_voting(logistic_model, gradient_boosting_model, x_test_scaled, y_test)

        logistic_predictions = logistic_model.predict_proba(x_combined_scaled)
        log_true_prob = [sublist[1] for sublist in logistic_predictions]
        combined_data["logistic_prediction_prob"] = log_true_prob

        gradient_boosting_predictions = gradient_boosting_model.predict_proba(x_combined_scaled)
        gradient_true_prob = [sublist[1] for sublist in gradient_boosting_predictions]
        combined_data["gradient_boosting_prediction_prob"] = gradient_true_prob
        combined_data["prediction"] = self.voting(
            logistic_model,
            gradient_boosting_model,
            x_combined_scaled
        )

        # read in user input for checked data
        #checked_data["prediction"] = checked_data["validated"]
        #combined_data['prediction'] =
        #combined_data['logistic_prediction'] | combined_data['gradient_boosting_prediction']
        # join sort and save
        #combined_data = pd.concat([checked_data, other_data], ignore_index=True)
        combined_data = combined_data.sort_values(by="IndexColumn")
        combined_data.to_csv(input_file, index=False)

        #print("\nLogistic Regression Coefficients:")
        #self.print_logistic_regression_coefficients(logistic_model, features)


    def train_logistic_regression(self, x_train, y_train):
        model = LogisticRegression()
        c_values = [0.01, 0.1, 1, 10, 100]
        best_recall = 0
        best_c = None
        for c in c_values:
            model = LogisticRegression(C=c, solver='liblinear')
            scores = cross_val_score(
                model,
                x_train,
                y_train,
                cv=5,
                scoring="recall"
            )
            avg_recall = scores.mean()
            if avg_recall > best_recall:
                best_recall = avg_recall
                best_c = c
                model.fit(x_train, y_train)
        return LogisticRegression(C=best_c, solver='liblinear').fit(x_train, y_train)


    def train_gradient_boosting(self, x_train, y_train):
        with parallel_backend('threading', n_jobs=1):
            param_grid = {
                'n_estimators': [50, 100, 150],
                'learning_rate': [0.01, 0.1, 0.2],
                'max_depth': [3, 5, 7],
                'min_samples_split': [2, 5, 10],
                'min_samples_leaf': [1, 2, 4]
            }
            scoring = make_scorer(recall_score)
            grid_search = GridSearchCV(
                estimator=GradientBoostingClassifier(),
                param_grid=param_grid,
                cv=5,
                n_jobs=-1,
                scoring=scoring
            )
            grid_search.fit(x_train, y_train)
        return grid_search.best_estimator_


    def evaluate_model(self, model, x_test, y_test):
        y_pred = model.predict(x_test)
        print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
        print("\nClassification Report:\n", classification_report(y_test, y_pred))
        accuracy = model.score(x_test, y_test)
        print("\nTest Set Accuracy:\n", accuracy)


    def print_logistic_regression_coefficients(self, model, features):
        coefficients = model.coef_[0]
        for feature, coefficient in zip(features, coefficients):
            print(f"{feature}: {coefficient}")


    def voting(self, logistic_model, gradient_boosting_model, x):
        logistic_proba = logistic_model.predict_proba(x)
        gradient_boosting_proba = gradient_boosting_model.predict_proba(x)

        vote_list = []
        for i in range(len(logistic_proba)):
            log_prob = logistic_proba[i][1]
            grad_prob = gradient_boosting_proba[i][1]
            if (log_prob > 0.5) or (grad_prob > 0.5):
                vote_list.append(True)
            else:
                vote_list.append(False)

        return vote_list

    def evaluate_voting(self, logistic_model, gradient_boosting_model, x_test, y_test):
        combined_predictions = self.voting(logistic_model, gradient_boosting_model, x_test)

        # evaluate the combined predictions
        print("\nModel Evaluation:")
        print("\nClassification Report:\n", classification_report(y_test, combined_predictions))
        tp = confusion_matrix(y_test, combined_predictions)[1][1]
        fn = confusion_matrix(y_test, combined_predictions)[1][0]

        recall = tp / (tp + fn)

        if recall <= 0.85:
            print("!WARNING!")
            print(f"Recall ({recall}) is less than 0.85")
            print("A signficant number of real peaks may be filtered from the data!")
            print("Please reattempt classification or use manual filtering if this is an issue")
            print("!WARNING!")
