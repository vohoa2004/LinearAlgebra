from matrix import Matrix


def input_matrix():
    rows = int(input("Nhap so hang: "))
    cols = int(input("Nhap so cot: "))

    matrix = []
    print(f"Nhap ma tran (moi dong cach nhau boi dau xuong hang):")
    for i in range(rows):
        row = list(map(int, input().split()))
        if len(row) != cols:
            raise ValueError(f"Dong {i + 1} phai co dung {cols} phan tu.")
        matrix.append(row)
    return matrix


print("Nhap ma tran A:")
matrixA = input_matrix()
A = Matrix(matrixA)

print("Ma tran A:")
print(A)

print("Nhap ma tran B:")
matrixB = input_matrix()
B = Matrix(matrixB)

print("Ma tran B:")
print(B)

F = A.transpose()
print("Ma tran chuyen vi cua ma tran A:")
print(F)

C = A.add(B)
print("A + B:")
print(C)

D = A.subtract(B)
print("A - B:")
print(D)

E = A.multiply(B)
print("A * B:")
print(E)

print("Ma tran don vi kich thuoc n x n")
n = int(input("Nhap n: "))
print(A.identity_matrix(n))

det_A = A.determinant()
print("Dinh thuc cua A:")
print(det_A)

inv_A = A.inverse()
print("Nghich dao cua A:")
print(inv_A)

rank_A = A.rank_matrix()
print("Hang cua A: ")
print(rank_A)

eigenvalues = A.eigenvalues()
print ("Cac tri rieng cua A: ")
print (eigenvalues)

eigenvectors = A.eigenvectors(eigenvalues)
print ("Cac vector rieng cua A: ")
print (eigenvectors)

P, D = A.diagonalize(eigenvalues, eigenvectors)
print("Ma tran duong cheo D sau khi cheo hoa ma tran A:\n", D)
print("D di cung voi P:\n", P)

isPositiveDefined = A.check_positive_define(eigenvalues)
print("Ma tran A la ma tran xac dinh duong" if isPositiveDefined else "Ma tran A khong phai la ma tran xac dinh duong")

vector_norm_l1 = Matrix.vector_norm_l1(A.matrix[0])
print("Chuan l1 cua vector A[0]: ", vector_norm_l1)

vector_norm_l2 = Matrix.vector_norm_l2(A.matrix[0])
print("Chuan l2 cua vector A[0]: ", vector_norm_l2)

frobenius_norm = A.matrix_norm_frobenius()
print("Chuan Frobenius cua Ma tran A: ", frobenius_norm)

# Testcase:
# A:
# 1 0 0
# -2 3 0
# 0 0 4



