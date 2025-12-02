
# ⚡ Charging Station Simulation

A complete simulation environment for studying charging stations based on **half-wave**, **full-wave**, and **center-tapped rectifiers**.  
The repository includes **MATLAB scripts**, **Simulink models**, and **Python implementations**, allowing you to analyze and compare rectifier behavior across multiple tools.

---

## 🚀 Features

### 🔧 MATLAB Features
- Simulation of:
  - Half-wave rectifier  
  - Full-wave rectifier  
  - Center-tapped full-wave rectifier  
- MATLAB scripts for testing and visualization  
- **Live script (`main.mlx`)**
- **Simulink implementations** of all rectifiers  
- **`simulink_params.m`** – parameter file used across all Simulink models  
- **`sweeper.m`** – automatically:
  - Sweeps through firing angles  
  - Computes average output voltage for all rectifiers  
  - Generates comparison plots  
  - Prepares all parameters in the workspace for Simulink models  

> **Run `sweeper.m` first**, then open or simulate any Simulink rectifier model.

### 🐍 Python Features
- Python implementations of:
  - Half-wave rectifier  
  - Full-wave center-tapped rectifier  
- Easy to modify and extend  
- Uses `numpy`, `matplotlib`, and `scipy`  

### General Features
- Multiple versions (MATLAB + Python) for flexibility  
- Easy to clone, run, and experiment  
- Good for learning, prototyping, and comparing methodologies  

---

## 📂 Repository Structure

```

/matlab
├── half_wave.slx
├── full_wave.slx
├── full_wave_center_tapped.slx
├── simulink_params.m
├── sweeper.m
├── main.mlx
└── (other MATLAB scripts)

python/
├── full_wave_rectifier_C_T.py
├── half_wave_rectifier.py

README.md

```


## 🛠️ Getting Started

### Requirements

#### MATLAB
- MATLAB R2021+ recommended  
- Simulink (required for .slx models)

#### Python
- Python 3.7+
- Required libraries:
  - `numpy`
  - `matplotlib`
  - `scipy`

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

4. Run any Simulink model:

   * `half_wave.slx`
   * `full_wave.slx`
   * `full_wave_center_tapped.slx`

5. Or open the Live Script:

   ```matlab
   main.mlx
   ```

---

## ▶️ Running the Python Version

```bash
cd python
python full_wave_rectifier_C_T.py
python half_wave_rectifier.py
```

### Optional: Create a virtual environment

```bash
python -m venv venv
source venv/bin/activate        # macOS / Linux
# or
venv\Scripts\activate           # Windows
```

Install dependencies:

```bash
pip install numpy matplotlib scipy
```

---

## 📊 MATLAB Example: Sweeper

Running `sweeper.m`:

* Iterates firing angle values
* Computes average output voltage for:

  * Half-wave
  * Full-wave
  * Center-tapped full-wave
* Generates comparison plots
* Makes parameters available for Simulink simulations

---

## 🎯 Why Both MATLAB and Python?

Having both environments allows you to:

* Validate analytical formulas using simulation
* Compare block-diagram Simulink results with Python numerical models
* Use MATLAB for visualization and model-based design
* Use Python for fast scripting, automation, and external integration

---

## 🤝 Contributing

Contributions are welcome!

1. Fork the repository
2. Create a new branch
3. Add your improvements
4. Submit a Pull Request

---

Happy Simulating! ⚡🔌
