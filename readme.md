```markdown
# charging-station

A simple simulation / example of a charging station based on full-wave and half-wave rectifiers.  
This repository provides two versions of the charging station implementation: one in MATLAB, and one in Python — you can run whichever you prefer.

## ⚙️ Features

- Full-wave rectifier charging station simulation  
- Half-wave rectifier charging station simulation  
- Both MATLAB and Python versions available  
- Easy to clone, run and experiment  

## 📂 Repository structure

```

/matlab     — MATLAB version of the charging station code
/python     — Python version of the charging station code
README.md   — this documentation file

````

## 🛠️ Getting Started

### Prerequisites

- For MATLAB version: a working installation of MATLAB.  
- For Python version: Python (preferably 3.7+) and required dependencies (if any — see below).

### Clone the repository

```bash
git clone https://github.com/Tesla-station/charging-station.git
cd charging-station
````

### Run the MATLAB version

1. Open MATLAB.
2. In MATLAB, open the project or navigate to the `matlab` folder. ([mathworks.com][1])
3. Run the main script (for example `main.m`, or as documented in the `matlab/` folder).

### Run the Python version

1. Navigate to the `python` folder:

   ```bash
   cd python
   ```
2. (Optional) Create and activate a virtual environment.

   ```bash
   python -m venv venv
   source venv/bin/activate   # Linux/MacOS
   # or `venv\Scripts\activate` on Windows
   ```
3. Install dependencies

    numpy
    matplot
    scipy

   Otherwise, ensure you have the needed libraries installed.

## 🧪 Example Usage

### MATLAB

```matlab
% In MATLAB:
cd('<path-to-repo>/charging-station/matlab');
test.m;
main.mlx % the live script
```

### Python

```bash
cd <path-to-repo>/charging-station/python
python full_wave_rectifier_C_T.py
python half_wave_rectifier.py
```

You should see output representing the charging station simulation (e.g. rectifier behavior, charging characteristics).

## 🎯 Why two implementations?

Providing both MATLAB and Python versions allows you to:

* Work in whichever environment you are more comfortable with
* Compare behavior between MATLAB and Python implementations
* Use MATLAB version for quick prototyping and visualization
* Use Python version for scripting, automation, or integration with other tools

## 📖 Contributing

Contributions are welcome! If you want to add features, fix bugs, or improve documentation:

1. Fork the repository.
2. Make your changes in a separate branch.
3. Test your changes (in MATLAB or Python).
4. Submit a Pull Request describing your changes.


Happy charging! ⚡

