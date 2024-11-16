A one-step script to convert GTH pseudopotential to ".upf" format using [CP2K's ATOM module](https://manual.cp2k.org/trunk/CP2K_INPUT/ATOM.html).

### Usage

1. Prepare a json file "input.json" containing the following infomation: (`cp2k_path` should be the path containing "exe" and "data" directory)
    ```json
    {
        "cp2k_path": "/home/pkgs/cp2k", 
        "element": "H",
        "xc": "LDA",
        "valence": 1,
        "prefix": "H-GTH-LDA"
    }
    ```
2. Run the following command:
    ```bash
    python3 gth2upf.py input.json
    ```
    It generates ".upf" pseudopotential file with the intermediate CP2K input (".inp") and output (".out") files in current directory.

