import subprocess
import os
import json

# Run pip-licenses and get the output in JSON format
result = subprocess.run(
    ["pip-licenses", "--format=json", "--with-license-file", "--no-license-path"], 
    capture_output=True, 
    text=True
)

print(result)

# Parse the JSON output
licenses = json.loads(result.stdout)

# Create a separate file for each module's license
for license_info in licenses:
    print(license_info)
    module_name = license_info["Name"]
    license_text = license_info["LicenseText"]

    # Define the file name as [moduleName]_license.txt
    file_name = f"{module_name}_license.txt"
    
    # Write the license text to the file
    with open(file_name, "w", encoding="utf-8") as file:
        file.write(license_text)

    print(f"Written license for {module_name} to {file_name}")
