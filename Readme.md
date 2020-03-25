# SHAKE128

This is an implementation of SHAKE128,  according to the FIPS202 standard. This project was made for a computer science course in cryptography, and should not be considered secure nor robust.

## Installation

Run `make` command to build the files.

## Test and usage 

You can launch the program using the following command:
```bash
./shake128 n < input
```
where `n` is the number of output bytes you want, and `input` being the file you want to hash. The output will be given in standard hexadecimal form.
 

You can test the software using this command, and you should get the corresponding output:
```console
➜ ./shake128 32 < /dev/null
7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26
```
