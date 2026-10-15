"""Entry point for the TP. Runs one exercise at a time, by name:

    python3 main.py 6.1a-i
    python3 main.py 6.2b-ii
    python3 main.py          # lists the available exercises

The names match the numbering in the TP PDF, so the command you type, the
section in the report and the file on disk all agree.
"""

import sys

from exercises.ex61 import a as ex61a, b as ex61b
from exercises.ex62 import b as ex62b

EXERCISES = {
    # 6.1 Analisis de mapeo
    "6.1a-i":   ex61a.i,
    "6.1a-ii":  ex61a.ii,
    "6.1a-iii": ex61a.iii,
    "6.1a-iv":  ex61a.iv,
    "6.1b":     ex61b.run,

    # 6.2 Llamado de variantes
    "6.2b-i":   ex62b.i,
    "6.2b-ii":  ex62b.ii,
}


def usage():
    print("usage: python3 main.py <exercise>")
    print("available exercises:")
    for name in EXERCISES:
        print(f"  {name}")


def main():
    if len(sys.argv) != 2:
        usage()
        return
    name = sys.argv[1]
    if name not in EXERCISES:
        print(f"'{name}' is not a known exercise.\n")
        usage()
        return
    EXERCISES[name]()


if __name__ == "__main__":
    main()
