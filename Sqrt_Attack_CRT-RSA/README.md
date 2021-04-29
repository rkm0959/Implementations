# Sqrt Attack for CRT-RSA

This repository implements an attack that factorizes `n` in time complexity `O(sqrt(min(d_p, d_q)))` (ignoring logs).

This problem was also submitted to picoCTF 2021 by me, with the easier version requiring `O(min(d_p, d_q))` complexity.

Despite the attack was mentioned in the paper "Twenty Years of Attacks on the RSA Cryptosystem", this attack is surprisingly hard to find the reference for, as it is also mentioned in the mathoverflow discussion in https://mathoverflow.net/questions/120160/attack-on-crt-rsa. 
