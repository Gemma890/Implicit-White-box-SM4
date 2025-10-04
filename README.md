# Implicit-White-box-SM4
A white-box SM4 implementation based on implicit function principles.
This repository contains Python scripts to generate a white-box SM4 implementation based on implicit function principles.
 following the method described in 
## Requirements
- Python 3 (version >= 3.7)
- BoolCrypt (version >= 0.1.1)
- SageMath equipped with CryptoMiniSat
## Reference Details
1. utils/equivalence.py, utils/functionalequations.py, utils/modularaddition.py and utils/utilities.py are from <a href="https://github.com/ranea/BoolCrypt">BoolCrypt</a>. Their lisences are MIT.
2. Third-party code used in this project follows their original open-source licenses. Complete license information can be found in <a href="https://github.com/ranea/BoolCrypt?tab=MIT-1-ov-file#readme">LICENSES.md</a>
 ## Usage(Linux)
This project provides a complete implementation of a white-box SM4 algorithm using hidden functions, featuring component generation and encryption verification.
### 1. Generating the implicit round functions
Firstly，we need to execute the <mark style="background-color: black; color: black;">equivalent_SM4.py</mark> script to get the relevant data for a fixed-key, implicit-function-based equivalent SM4 algorithm. The transformations represented by this data possess a dual significance: they function as equivalent variants of the conventional SM4 algorithm and constitute the critical raw material required for the subsequent white-box implementation.
Then if you want to figure out whether the implicit-function-based equivalent SM4 algorithm we get just now runs in the same way as original SM4 algorithm,  you could execute <mark style="background-color: black; color: black;">equivalent_SM4_encrypt.py</mark> to test.
### 2. Generating encodings
To protect the intermediate data of the encryption process within the white-box context, we need to run <mark style="background-color: black; color: black;">generate_wbsm4_tool.py</mark> to generate various encodings. These encodings are then applied to the transformations from the previous equivalent algorithm, resulting in a white-box SM4 implementation based on implicit function structure. Due to the involvement of extensive large-scale matrix computations and pseudorandom number generation, this process is highly demanding in terms of both time and memory consumption.
### 3. Evaluating the implicit white-box implementation of SM4
The final outcome is a validated fixed-key white-box implementation by executing <mark style="background-color: black; color: black;">evaluate_wbsm4.py</mark> script, ready for direct use in encryption tasks, which confirms the correctness and efficiency of our code.
