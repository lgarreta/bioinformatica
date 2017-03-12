
# Program to compute the global and local pairwise alignment
# It implements the Smith-Waterman Algorithm for global
# alignments and The Smith-Waterman Algorithm for local 
# alignments.


def main ():
    seq1 = "AGCGTAG"
    seq2 = "CTCGTC"
    #scoringPenalties = {"match" : 1, "mismatch" : 0, "gap" : -1}
    scoringPenalties = {"match" : 10, "mismatch" : -5, "gap" : -7}

    print "\n>>>Global Alignment of %s and %s" % (seq1, seq2)
    matrix, optimalScore, optimalLocalization = global_align (seq1,seq2, scoringPenalties, getMaxScoreGlobal)
    printMatrix (seq1, seq2, matrix)
    print "optimal score:", optimalScore
    print "optimal localization:", optimalLocalization

    print "\n>>>Local Alignment of %s and %s" % (seq1, seq2)

    matrix, optimalScore, optimalLocalization = local_align (seq1,seq2, scoringPenalties, getMaxScoreLocal)
    printMatrix (seq1, seq2, matrix)
    print "optimal score:", optimalScore
    print "optimal localization:", optimalLocalization

# Construct the alignment matrix */
def make_matrix(seq1, seq2):
	nColumns = len (seq1) + 1
	nFilas = len (seq2) + 1

	"""Creates a sizex by sizey matrix filled with zeros."""
	return [[0]*nColumns for i in xrange(nFilas)]

# Print the matrix with the partial scores
def printMatrix(seq1, seq2, matrix):
  seq1 = "-" + seq1
  seq2 = "-" + seq2

  print "-",
  for j in seq1:
	print j.rjust(3),
  print ("")
  
  for i in range (0, len (matrix)):
	print seq2 [i],
	for j in range (0, len (matrix[i])):
	  print str (matrix [i][j]).rjust(3),

	print ""

# Return the partial maximal score for a global alignment
def getMaxScoreGlobal (A, i, j, seq1, seq2, scoringPenalties):
    maxScore = max (A[i][j-1] + scoringPenalties ["gap"],
                 A[i-1][j] + scoringPenalties ["gap"],
                 A[i-1][j-1] + (scoringPenalties ["match"] if seq2[i-1]==seq1[j-1] else scoringPenalties ["mismatch"])
            )
    return maxScore

# Return the partial maximal score for a local alignment
def getMaxScoreLocal (A, i, j, seq1, seq2, scoringPenalties):
    maxScore = max (A[i][j-1] + scoringPenalties ["gap"],
                 A[i-1][j] + scoringPenalties ["gap"],
                 A[i-1][j-1] + (scoringPenalties ["match"] if seq2[i-1]==seq1[j-1] else scoringPenalties ["mismatch"]),
                 0
            )
    return maxScore
 
# Construct a global alignment
def global_align(seq1, seq2, scoringPenalties, scoringFunction):
	"""Do a local alignment between x and y"""
	# create a zero-filled matrix
	A = make_matrix (seq1, seq2)
	A[0][0] = 0
	  
	for j in range (1, len(seq1)+1):
	  A[0][j] += scoringPenalties ["gap"] + A[0][j-1]
	
	for i in range (1, len(seq2)+1):
	  A[i][0] += scoringPenalties["gap"] + A[i-1][0]

	return align (seq1, seq2, A, scoringPenalties, scoringFunction) 

def local_align(seq1, seq2, scoringPenalties, scoringFunction):
	"""Do a local alignment between x and y"""
	# create a zero-filled matrix
	A = make_matrix (seq1, seq2)
	A[0][0] = 0
	  
	for j in range (1, len(seq1)+1):
	  A[0][j] = 0
	
	for i in range (1, len(seq2)+1):
	  A[i][0] = 0

	return align (seq1, seq2, A, scoringPenalties, scoringFunction) 
 

# General method for alignments using as a parameter the scoring
# penalties and the scoring function
def align(seq1, seq2, A, scoringPenalties, scoringFunction):
    best = 0
    ptloc = (0,0)
    # fill in A in the right order
    for i in xrange(1, len(seq2)+1):
        for j in xrange(1, len(seq1)+1):
            # the local alignment recurrance rule:
            A[i][j] = scoringFunction (A, i, j, seq1, seq2,  scoringPenalties)
            # track the cell with the largest score
            if A[i][j] >= best:
               best = A[i][j]
               optloc = (i,j)
  # return the opt score and the best location
    return A,best, optloc

#--------------------------------------------------------------------
# Main
#--------------------------------------------------------------------
if __name__ == "__main__":
	main ()

