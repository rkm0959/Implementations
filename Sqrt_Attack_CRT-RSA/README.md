# Sqrt Attack for CRT-RSA

This repository implements an attack that factorizes `n` in time complexity `O(sqrt(min(d_p, d_q)))` (ignoring logs).

This problem was also submitted to picoCTF 2021 by me, with the easier version requiring `O(min(d_p, d_q))` complexity.

Despite the attack was mentioned in the paper "Twenty Years of Attacks on the RSA Cryptosystem", this attack is surprisingly hard to find the reference for, as it is also mentioned in the mathoverflow discussion in https://mathoverflow.net/questions/120160/attack-on-crt-rsa. 

The latter part of multipoint evaluation can be optimized with Chirp-Z transform. This is because the evaluation points are geometric. For further details, consult the writeup of ezRSA+++ at https://rkm0959.tistory.com/264. The implementation is also on https://github.com/rkm0959/Cryptography_Writeups/blob/main/2022/0CTF/ezRSA%2B%2B%2B.py
