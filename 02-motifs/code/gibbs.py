#!/usr/bin/python
from Bio import motifs
import random

def logodd (matPWM, back, S):
    scoreModelo = 1
    for k, i in enumerate (S):
        scoreModelo = scoreModelo * (matPWM [i][k] + 0.01)

    return scoreModelo

W=4
SCOREINICIAL = -1
ERROR = 0.5

seqsTmp = open ("sequences.txt").readlines ()
sequences = [x.strip() for x in seqsTmp]
T = len (sequences)
posicionesRandom = range (0,T)
random.shuffle (posicionesRandom)
print " > Posiciones aleatorias para escoger Sh:\n", posicionesRandom

N = len (sequences [0])

mt = motifs.create (sequences)

print "> Motivo:\n", mt

background = {"A":0.25,"C":0.25,"G":0.25,"T":0.25}

alineamientos = []
for seq in sequences:
    n = len (seq)
    i = random.randint (0, n-W)
    subseq = seq [i:i+W]
    alineamientos.append (subseq)

print "> Alineamientos:\n", alineamientos

mtaln = motifs.create (alineamientos)
print "> Motif alineamientos:\n", mtaln
pwm = mtaln.pwm

# Seleccion al azar Sh
iSh = posicionesRandom [0]
seqSh = sequences [iSh]

listShs = []
listShsScores = []
maxIndice = 0
maxScore = -1
maxSubseq = ""
for i in range (0, N-W+1):
    subseq = seqSh [i:i+W]
    listShs.append (subseq)
    score = logodd (pwm, background, subseq)
    if score > maxScore: 
        maxIndice = i
        maxScore = score
        maxSubseq = subseq

print "> Todos los motivios de Sh", seqSh, " :\n", listShs
print "> El mejor: %s, score: %s y pos: %s" % (maxSubseq, maxScore, maxIndice)
print pwm

diferenciaScores = abs (SCOREINICIAL - maxScore)
if (diferenciaScores < ERROR):
    print "El consenso del motivo es: %s", maxScore
    print "La PWM del motivo es: %s", pwm


