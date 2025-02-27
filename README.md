

# V4C - Virtual 4C Analysis for Hi-C Data
![BuildStatus](https://github.com/zoechen0717/BMI203-HW4-Clustering/workflows/badge.svg?event=push)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.6%2B-blue)](https://www.python.org/downloads/)

## 🚀 Overview
**V4C** is a Python package designed for **Hi-C Virtual 4C analysis**, allowing users to:
- **Extract Hi-C contact frequencies** from `.mcool` files at multiple resolutions.
- **Plot Virtual 4C contact profiles** with customizable parameters.
- **Compare multiple Hi-C datasets** to visualize differences across conditions.

The package supports:
✅ **Command-line usage** and **Python API**  
✅ **Custom genomic coordinates or default promoter regions (hg38/hg19)**  
✅ **Multiple `.mcool` inputs & resolutions**  
✅ **ICE balancing and Min-Max normalization**  
✅ **Custom upstream/downstream window size**  

---

## 📦 Installation

### 1️⃣ **Install via GitHub**
```bash
pip install git+https://github.com/zoechen0717/V4C.git
```

### 2️⃣ **Clone & Install Locally**
```bash
git clone https://github.com/zoechen0717/V4C.git
cd V4C
pip install -e .
```

### 3️⃣ **Check Installation**
```bash
v4c-extract --help
```

---

## 🛠️ Dependencies
`V4C` requires the following Python libraries:
```txt
numpy
pandas
matplotlib
cooler
argparse
```
Install them manually if needed:
```bash
pip install numpy pandas matplotlib cooler argparse
```

---

# 👉 Command-Line Usage

### **1️⃣ Extract Hi-C Contact Frequencies**
Extract contact frequencies from `.mcool` files at **multiple resolutions**.

```bash
v4c-extract --mcool sample1.mcool sample2.mcool \
            --resolution 5000 10000 \
            --coords chr17:45878152-46000000 \
            --genome hg38 \
            --bed custom_regions.bed \
            --flank 50000 \
            --balance True \
            --scale True \
            --output extracted_data.tsv
```

#### **✅ Explanation of Options**
| Argument | Description |
|----------|-------------|
| `--mcool` | One or multiple `.mcool` files to process |
| `--resolution` | Resolution(s) to extract (e.g., 5000, 10000) |
| `--coords` | Genomic coordinates (e.g., `chr17:45878152-46000000`) |
| `--genome` | `hg38` or `hg19` (default promoter locations) |
| `--bed` | Custom BED file with genomic coordinates (optional) |
| `--flank` | Flanking region in bp (default: `50000`) |
| `--balance` | Apply ICE normalization (`True/False`) |
| `--scale` | Min-Max scale between 0-1 (`True/False`) |
| `--output` | Output file name (default: `extracted_data.tsv`) |

---

### **2️⃣ Plot Virtual 4C Data**
Visualize extracted **Hi-C contact frequency profiles**.

```bash
v4c-plot --input extracted_data.tsv --ylim 0.4 --flank 50000
```

#### **✅ Explanation of Options**
| Argument | Description |
|----------|-------------|
| `--input` | Extracted `.tsv` file from `v4c-extract` |
| `--ylim` | Maximum y-axis value (default: `0.4`) |
| `--flank` | Flanking region in bp |

---

### **3️⃣ Compare Multiple Hi-C Datasets**
Compare **multiple Hi-C datasets** in the same figure.

```bash
v4c-compare --inputs extracted_data1.tsv extracted_data2.tsv \
            --ylim 0.4 \
            --scale True
```

#### **✅ Explanation of Options**
| Argument | Description |
|----------|-------------|
| `--inputs` | One or multiple `.tsv` files from `v4c-extract` |
| `--ylim` | Maximum y-axis value |
| `--scale` | Min-Max scale between 0-1 (`True/False`) |

---

# 📂 Example Output

### **Extracted Data (`extracted_data.tsv`)**
```
mcool    res    chrom    start      end        contact_freqs
sample1  5000   chr17    45878152   46000000   0.12,0.24,0.34,0.15,0.08
sample2  5000   chr17    45878152   46000000   0.10,0.22,0.32,0.14,0.07
```

### **Virtual 4C Plot**
![V4C Plot Example](https://upload.wikimedia.org/wikipedia/commons/a/a4/4C-seq_diagram.png)

### **Comparison Plot**
![Comparison Example](https://upload.wikimedia.org/wikipedia/commons/thumb/3/3e/DNA_replication_en.svg/1200px-DNA_replication_en.svg.png)

---

## 🚀 Development & Contribution

### **1️⃣ Clone the Repository**
```bash
git clone https://github.com/yourusername/V4C.git
cd V4C
pip install -e .
```

### **2️⃣ Run Tests**
```bash
pytest tests/
```

### **3️⃣ Submit a Pull Request**
- Fork the repository
- Create a new branch (`feature-x`)
- Make your changes and push to GitHub
- Open a Pull Request

---

## 💜 License
**V4C** is released under the MIT License. See [LICENSE](LICENSE) for details.

---

## 📞 Contact
For issues and discussions, visit the [GitHub Issues](https://github.com/yourusername/V4C/issues) page.

📧 **Email**: your.email@example.com  
🏡 **GitHub**: [yourusername](https://github.com/yourusername)
