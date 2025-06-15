#testing
import sys
import os
from interface import *
from sequence import *
from lcsfinder import *
from sequencealignment import *
from sequencedatabase import *

### Main function
def main():
    """Main function to run the Gene LCS application.
    :return: int - exit code
    """
    ## ve se o ficheiro de entrada foi fornecido como argumento e se existe
    if len(sys.argv) >= 2 and os.path.exists(sys.argv[1]):
        app = GeneLCSApp(sys.argv[1])
        app.mainloop()
    else:
        print("Error: Input file not provided or non existent\nUsage: python .\\projeto_1 data\\<INPUT_FILE_PATH>")
        return -1

# Redirects to main
if __name__ == "__main__":
    main()