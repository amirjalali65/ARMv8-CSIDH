# ARMv8-CSIDH
Optimized and Constant-time implementation of Post-Quantum Commutative Supersingular Isogeny Diffie-Hellman (CSIDH) key exchange on ARMv8 Processors.

This repository contains an highly-optimized implementation of CSIDH on ARMv8 processors. The finite field arithmetic is designed and engineered for the p511 prime proposed in the original [CSIDH](https://eprint.iacr.org/2018/383.pdf) scheme by Castryck et al. 

This library implements CSIDH in both constant-time and variable-time. The constant-time implementation is designed and developed by Amir Jalali. The variable-time implementation is developed based on the proof-of-concept implementation of CSIDH authors. 

## Building Binaries
### Cross Compliation for ARMv8 on Ubuntu
ARMv8 executables can be generated using cross-compilation on Linux. There are different methods for cross-compilation. An easy approach is to install `gcc-aarch64-linux-gnu` package by executing:
```sh
$  sudo apt-get install gcc-aarch64-linux-gnu
```
After installation, use the following command to generate the ARMv8 executables:
### Variable-time Executable
Simply use `make` in the terminal:
```sh
$ make 
```
### Constant-time Executable
Use the following command in the terminal:
```sh
$ make CONSTANT=TRUE 
```
### Constant-time with uniform variable-time ladder
In order to improve the constant-time CSIDH performance, we can replace the constant-time Montgomery ladder inside the scheme with the uniform variable-time ladder. This rises some security concerns regarding side-channel attacks, however the performance improvement is significant. 
To generate a constant-time with uniform (variable-time) ladder, use the following command in the terminal:
```sh
$ make CONSTANT=TRUE FASTLADDER=TRUE
```


The generated executable is `CSIDH_TEST` and can be run on ARMv8 cores.


## Contributors
The constant-time implementation as well as optimized finite field arithmetic are designed and developed by Amir Jalali (ajalali2016@fau.edu).
The variable-time implementation and key validation is designed with minor modifications based on the CSIDH proof-of-concept implementation by Castryck et al. The field arithmetic implementation is designed for ARMv8 processors.





