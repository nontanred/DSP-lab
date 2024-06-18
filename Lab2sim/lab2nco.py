
# NCO function for 3rd part of experiment #2
# input:
#       n - discrete time argument [int]
#       k - step size as frequency multiplier [int]
#       sin_array - array of 64 samples of sine wave period [int as Q15]
# output:
#       sin - sin(2pi * n * k) [int as Q15]
def NCO(n, k, sin_array):
    sin_array_idx = (n) % 64
    sin = sin_array[sin_array_idx]
    return sin

# sin_array - array of 64 samples of sine wave period [int as Q15]
def sin_array():
    sin_array64 = [0,
                 3212,
                 6393,
                 9512,
                 12539,
                 15446,
                 18204,
                 20787,
                 23170,
                 25329,
                 27245,
                 28898,
                 30273,
                 31356,
                 32137,
                 32609,
                 32767,
                 32609,
                 32137,
                 31356,
                 30273,
                 28898,
                 27245,
                 25329,
                 23170,
                 20787,
                 18204,
                 15446,
                 12539,
                 9512,
                 6393,
                 3212,
                 0,
                 -3212,
                 -6393,
                 -9512,
                 -12539,
                 -15446,
                 -18204,
                 -20787,
                 -23170,
                 -25329,
                 -27245,
                 -28898,
                 -30273,
                 -31356,
                 -32137,
                 -32609,
                 -32767,
                 -32609,
                 -32137,
                 -31356,
                 -30273,
                 -28898,
                 -27245,
                 -25329,
                 -23170,
                 -20787,
                 -18204,
                 -15446,
                 -12539,
                 -9512,
                 -6393,
                 -3212]
    return sin_array64