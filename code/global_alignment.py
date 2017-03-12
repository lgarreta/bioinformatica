
# Program to compute the global and local pairwise alignment
# It implements the Smith-Waterman Algorithm for global
# alignments 

# Main function where the program starts
def main ():
    seq1 = "ACGT"
    seq2 = "ACCGCT"
    scoringPenalties = {"match" : 2, "mismatch" : -1, "gap" : -2}
    #scoringPenalties = {"match" : 10, "mismatch" : -5, "gap" : -7}
    print "\n>>> Global Alignment of %s and %s" % (seq1, seq2)
    matrix = global_align (seq1,seq2, scoringPenalties)
    printMatrix (seq1, seq2, matrix)

# Construct the alignment matrix filled with zeros 
def make_matrix(seq1, seq2):
    nColumns = len (seq1) + 1
    nFilas = len (seq2) + 1
    return [[0]*nColumns for i in xrange(nFilas)]

# Construct a global alignment
def global_align(seq1, seq2, scoringPenalties):
    # create a zero-filled matrix
    A = make_matrix (seq1, seq2)
    A[0][0] = 0

    # Initialize the first row with gap penalties 
    for i in range (1, len(seq2)+1):
      A[i][0] += scoringPenalties["gap"] + A[i-1][0]
    # Initialize the first column with gap penalties  
    for j in range (1, len(seq1)+1):
      A[0][j] += scoringPenalties ["gap"] + A[0][j-1]
    
    # fill in A in the right order
    for i in xrange(1, len(seq2)+1):
	for j in xrange(1, len(seq1)+1):
	    typeSubstitution = "match" if seq2[i-1] == seq1[j-1] else "mismatch"

	    diagonalScore   = A[i-1][j-1] + scoringPenalties [typeSubstitution]
	    horizontalScore = A[i][j-1] + scoringPenalties ["gap"]
	    verticalScore   = A[i-1][j] + scoringPenalties ["gap"]

	    A[i][j] = max (horizontalScore, verticalScore, diagonalScore)

    return A

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


#--------------------------------------------------------------------
# Main
#--------------------------------------------------------------------
if __name__ == "__main__":
	main ()

