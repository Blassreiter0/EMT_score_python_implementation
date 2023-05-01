# EMT_score_python_implementation

This project takes a GSEID as input and downloads the GSE dataset and its annotation table. It then implements certain algorithms to generate 76GS score, KS Score, and MLR Score for each GSM ID. Then it generates a plot of which correlates the scores obtained from the 3 methods.


## Installation

Make sure you have matlab runtime R2023a installed on your device. It can be installed from:-

https://in.mathworks.com/products/compiler/matlab-runtime.html


Clone the repsitory
```
git clone https://github.com/Blassreiter0/EMT_score_python_implementation.git
```

Create a virtual environment:
```
python -m venv venv
```

Activate the virtual environment:
```
source venv/bin/activate
```

Install the requirements:
```
pip install -r requirements.txt
```


## Usage

go inside the repository using the command:
cd .\EMT_score_python_implementation\

Open the main.py file from any IDE. Run it and Then enter any GSEID(for Eg. GSE21422).


## Credits

This code is the python implementation of work from github projects:-<br>
https://github.com/priyanka8993/EMT_score_calculation<br>
https://github.com/sushimndl/EMT_Scoring_RNASeq<br>
<br>
for more information:<br><br>
Chakraborty P, George JT, Tripathi S, Levine H and Jolly MK (2020) Comparative Study of Transcriptomics-Based Scoring Metrics for the Epithelial-Hybrid-Mesenchymal Spectrum. Front. Bioeng. Biotechnol. 8:220. doi: 10.3389/fbioe.2020.00220
