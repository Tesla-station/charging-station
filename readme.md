
# ⚡ Charging Station Simulation

A simulation toolbox for studying rectifier-based charging stations.  
This repository now includes MATLAB scripts, Simulink models, and Python implementations for analyzing half-wave, full-wave, and center-tapped rectifiers.

---

## ✨ New Features (Updated)

The MATLAB section now includes:

- **Simulink models** for:
  - Half-wave rectifier  
  - Full-wave rectifier  
  - Center-tapped full-wave rectifier  

- **`simulink_params.m`**  
  Centralized parameter file for all Simulink models.

- **`sweeper.m`**  
  A MATLAB script that:
  - Sweeps firing angles  
  - Computes the **average output voltage** for each rectifier  
  - Generates comparative plots  

> **Important:**  
> Run **`sweeper.m` first**, because it loads the required parameters into the MATLAB workspace.  
> After that, you may run any Simulink rectifier model individually.



## 📁 Repository Structure


/matlab
├── half_wave.slx
├── full_wave.slx
├── full_wace_center_tapped.slx
├── simulink_params.m
├── sweeper.m
└── (other MATLAB scripts)

python/
README.md


## 🛠️ Getting Started

### Requirements

- **MATLAB** (with Simulink)
- **Python 3.7+** (for the Python implementation)
- Python libraries:
  - numpy  
  - matplotlib  
  - scipy  

---

## ▶️ Running the MATLAB Version

1. Open MATLAB.
2. Navigate to the MATLAB folder:

   ```matlab
   cd('<path-to-repo>/matlab')


3. **Run the parameter and sweep script first:**

   ```matlab
   sweeper
   ```

   This populates the workspace with parameters used by all Simulink models.

4. To run a specific rectifier model, simply open and run it:

   * `half_wave_rectifier.slx`
   * `full_wave_rectifier.slx`
   * `center_tapped_rectifier.slx`

5. Or run any additional MATLAB scripts provided in the folder.

---

## ▶️ Running the Python Version

```bash
cd python
python full_wave_rectifier_C_T.py
python half_wave_rectifier.py
```

Optional: use a virtual environment.

```bash
python -m venv venv
source venv/bin/activate     # Linux/macOS
# or venv\Scripts\activate   # Windows
```

Install dependencies:

```bash
pip install numpy matplotlib scipy
```

---

## 📊 MATLAB Example: Sweeper Plot

`sweeper.m` automatically:

* Iterates over firing angle values
* Computes the average output voltage for:

  * Half-wave
  * Full-wave
  * Center-tapped
* Displays a comparison plot
* Stores results in workspace variables for reuse in Simulink

---

## 🎯 Purpose of This Project

Including both MATLAB/Simulink and Python implementations allows you to:

* Compare analytic vs. simulation-based results
* Use Simulink for block-level modeling
* Use MATLAB scripts for parameter sweeps and automated plots
* Use Python for scripting, automation, or external tooling

---

## 🤝 Contributing

Contributions are welcome!

1. Fork the repository
2. Create a new branch
3. Add your improvements
4. Open a Pull Request

---

Happy Simulating! ⚡🔌