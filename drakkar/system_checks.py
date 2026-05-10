import os
import subprocess
import sys

from drakkar.output import Text as RichText, print, prompt

def is_snakemake_locked(workdir: str) -> bool:
    locks_dir = os.path.join(workdir, ".snakemake", "locks")
    return os.path.isdir(locks_dir) and len(os.listdir(locks_dir)) > 0

def check_screen_session():
    """Checks if the script is running inside a 'screen' session. If not, warns the user and asks for confirmation."""
    if "STY" not in os.environ:
        print("\n ⚠️   WARNING: You are not running this script inside a 'screen' session.")
        print("     Running long processes outside of 'screen' may cause issues if your session is disconnected.")
        if RichText is None:
            print("\n 📌   To start a screen session, use: `screen -S mysession`")
        else:
            screen_message = RichText("\n 📌   To start a screen session, use: ")
            screen_message.append("`screen -S mysession`", style="drakkar.code")
            print(screen_message)
        print("         Then run this script inside the screen session.\n")

        # Prompt user to continue
        while True:
            user_input = prompt("    👉 Type '1' to ignore this warning and continue, or Ctrl+C to exit: ")
            if user_input.strip() == "1":
                break  # Continue execution
            else:
                print("     Invalid input. Please type '1' to continue.")
