from rich.console import Console

console = Console(highlight=False)

def info(msg):
    console.print(f"{msg}")

def warn(msg):
    console.print(f"[yellow]⚠  ATTENTION \n{msg}[/yellow]")

def error(msg):
    console.print(f"[red]❌ [b]ERROR[/b]  \n{msg}[/red]")

def success(msg):
    console.print(f"[green]✅ [b]STEP COMPLETE[/b] \n{msg}[/green]")
