from pprint import pprint
import numpy as np

LEFT = "--"
UP = "||"
DIAG = "\\"


def listtostring(list):
    string = ""
    for e in list:
        string += e

    return string


def dir(left, diag, up, x, y):
    if left > diag and left > up:
        return left, LEFT
    if up > diag and up > left:
        return up, UP
    if diag > left and diag > up:
        return diag, DIAG
    if left == up and up > diag:
        if x < y:
            return up, UP
        else:
            return left, LEFT

    if left == diag and diag > up:
        if x > y:
            return left, LEFT
        else:
            return diag, DIAG

    if diag == up and up > left:
        if x < y:
            return up, UP
        else:
            return diag, DIAG

    if left == up == diag:
        if x > y:
            return left, LEFT
        elif x < y:
            return up, UP
        else:
            return diag, DIAG

    return diag, DIAG


def needleman(sequ2, sequ1, match=3, mismatch=-2, gap=-1):
    score = np.full((len(sequ1) + 1, len(sequ2) + 1), 0)
    direction = np.full((len(sequ1) + 1, len(sequ2) + 1), "NONE")

    for y in range(0, len(score)):
        for x in range(0, len(score[0])):
            if y == 0:
                direction[y][x] = LEFT
                score[y][x] = -1 * x
            if x == 0:
                direction[y][x] = UP
                score[y][x] = -1 * y

    direction[0][0] = DIAG

    for y in range(1, len(score)):
        for x in range(1, len(score[0])):
            left = score[y][x - 1] + gap
            up = score[y - 1][x] + gap
            diag = 0
            if sequ1[y - 1] == sequ2[x - 1]:
                diag = score[y - 1][x - 1] + match
            else:
                diag = score[y - 1][x - 1] + mismatch

            temp = dir(left, diag, up, x, y)
            score[y][x] = temp[0]
            direction[y][x] = temp[1]
    return score, direction


def needleman_backtrack(directions, sequ2, sequ1):
    sequ2 = "-" + sequ2
    sequ1 = "-" + sequ1
    x = len(directions) - 1
    y = len(directions[0]) - 1
    gen_sequenz = [[], []]
    while x > 0 or y > 0:
        if directions[x][y] == DIAG:
            gen_sequenz[0].append(sequ1[x])
            gen_sequenz[1].append(sequ2[y])
            x -= 1
            y -= 1
        elif directions[x][y] == UP:
            gen_sequenz[0].append(sequ1[x])
            gen_sequenz[1].append("-")
            x -= 1
        elif directions[x][y] == LEFT:
            gen_sequenz[0].append("-")
            gen_sequenz[1].append(sequ2[y])
            y -= 1

    seq1 = gen_sequenz[1]
    seq2 = gen_sequenz[0]
    seq1 = seq1[::-1]
    seq2 = seq2[::-1]
    seq1 = listtostring(seq1)
    seq2 = listtostring(seq2)
    return seq1, seq2


sequ1 = "ACCGAAGTAC"
sequ2 = "AGAGGTAC"
score, direction = needleman(sequ1, sequ2, match=3, mismatch=-2, gap=-3)
pprint(score, indent=1, width=300)
pprint(direction, indent=1, width=300)
print(needleman_backtrack(direction, sequ1, sequ2), "\n")

score, direction = needleman("ACCGAAGTAC", "AGAGGTAC", match=1, mismatch=-4, gap=-1)
pprint(score, indent=1, width=300)
pprint(direction, indent=1, width=300)
print(needleman_backtrack(direction, sequ1, sequ2), "\n")

score, direction = needleman("ACCGAAGTAC", "AGAGGTAC", match=2, mismatch=-5, gap=-3)
pprint(score, indent=1, width=300)
pprint(direction, indent=1, width=300)
print(needleman_backtrack(direction, sequ1, sequ2), "\n")
