import math

class Matrix:
    def __init__(self, matrix=None):
        if matrix is not None:
            self.matrix = matrix
            self.R = len(matrix)
            self.C = len(matrix[0])
        else:
            self.matrix = []
            self.R = 0
            self.C = 0

    # - Chuyển vị ma trận
    def transpose(self):
        result = [
                    [self.matrix[i][j]
                    for i in range(self.R)] # so hang
                    for j in range(self.C) # so cot
                  ] # hang thanh cot, cot thanh hang
        return Matrix(result)

    # - Phép cộng, trừ, nhân 2 hoặc nhiều ma trận
    def add(self, other):
        if self.shape() != other.shape():
            raise ValueError("Cac ma tran phai co cung kich thuoc moi cong duoc")
        result = [[self.matrix[i][j] + other.matrix[i][j]
                    for j in range(self.C)]
                    for i in range(self.R)]
        return Matrix(result)

    def subtract(self, other):
        if self.shape() != other.shape():
            raise ValueError("Cac ma tran phai co cung kich thuoc moi tru duoc")
        result = [[self.matrix[i][j] - other.matrix[i][j]
                   for j in range(self.C)]
                  for i in range(self.R)]
        return Matrix(result)

    def multiply(self, other):
        if self.C != other.R: # so cot cua ma tran A != so dong ma tran B => ko nhan dc
            raise ValueError(
                "Ma tran thu nhat phai co so cot bang so hang cua ma tran thu hai moi nhan duoc"
            )
        result = [
                    [sum(self.matrix[i][k] * other.matrix[k][j]
                        for k in range(other.R)) # so dong cua ma tran B
                        for j in range(other.C) # so cot cua ma tran B
                    ]
                        for i in range(self.R)
                ] # so dong cua ma tran A
        return Matrix(result)

    # - Lấy một ma trận đơn vị n x n voi n la so hang cua self matrix
    @staticmethod
    def identity_matrix(n):
        result = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        return Matrix(result)

    # - Tính định thức một ma trận
    # Cong thuc: det A = (-1)^(i + j) * a_ij * det(A_ij)
    def determinant(self):
        if self.R != self.C:
            raise ValueError("Ma tran phai vuong moi tinh dinh thuc duoc")
        return self._determinant_recursive()

    def _determinant_recursive(self):
        if self.C == 1:
            return self.matrix[0][0]
        if self.C == 2:
            return self.matrix[0][0] * self.matrix[1][1] - self.matrix[0][1] * self.matrix[1][0] # det cuama tran 2 x 2

        # khai trien theo dong dau tien => i + j = j => j di tu 0 toi het so rows cua ma tran
        det = 0
        for c in range(len(self.matrix)): # square matrix: R = C
            sub_matrix = Matrix([row[:c] + row[c + 1:] for row in self.matrix[1:]]) # phan bu dai so, loai bo cot thu c, hang 0
            det += ((-1) ** c) * self.matrix[0][c] * sub_matrix._determinant_recursive() # (-1)^(i+j) * a_ij * det(A_ij) voi i = 0, j = c
        return det

    # - Nghịch đảo một ma trận
    # Cong thuc: inverse A = (adj A)/(det A)
    def inverse(self):
        det = self.determinant()
        if det == 0:
            raise ValueError("Dinh thuc bang 0 nen ma tran khong kha nghich")
        return Matrix(self._adjugate(self.matrix, det))

    def _adjugate(self, matrix, det):
        cofactors = []
        for r in range(self.R):
            cofactor_row = []
            for c in range(self.R):
                minor = Matrix([row[:c] + row[c + 1:] for row in (matrix[:r] + matrix[r + 1:])])
                cofactor_row.append(((-1) ** (r + c)) * minor._determinant_recursive())
            cofactors.append(cofactor_row)
        cofactors = Matrix(cofactors).transpose().matrix
        return [[cofactors[r][c] / det for c in range(len(cofactors))] for r in range(len(cofactors))]

    # - Tính Rank của một ma trận
    def rank_matrix(self):
        matrix = self.matrix
        rank = len(matrix[0])
        for row in range(0, rank, 1):
            if matrix[row][row] != 0:
                for col in range(0, len(matrix), 1):
                    if col != row:

                        multiplier = (matrix[col][row] /
                                      matrix[row][row])
                        for i in range(rank):
                            matrix[col][i] -= (multiplier *
                                               matrix[row][i])
            else:
                reduce = True

                for i in range(row + 1, len(matrix), 1):
                    if matrix[i][row] != 0:
                        self._swap(row, i, rank)
                        reduce = False
                        break

                if reduce:
                    rank -= 1
                    for i in range(0, len(matrix), 1):
                        matrix[i][row] = matrix[i][rank]
                row -= 1
        return rank

    def _swap(self, row1, row2, col):
        for j in range (col):
            temp = self.matrix[row1][j]
            self.matrix[row1][j] = self.matrix[row2][j]
            self.matrix[row2][j] = temp

    # - Tính trị riêng của một ma trận

    # Cach 1:
    # Da thuc dac trung: P_An(eigen_value) = det(eigen_value * I_n - A_n)
    # eigen value = n0 cua Da thuc dac trung

    def eigenvalues(self):
        if self.C != self.R:
            raise ValueError("Ma tran phai vuong moi tim duoc tri rieng")
        eigenvalues = []
        for lambd in range(-100, 100):
            lambda_matrix =  self.create_lambda_matrix(lambd)
            if abs(lambda_matrix._determinant_recursive()) < 1e-6:
                eigenvalues.append(lambd)
        return eigenvalues

    def create_lambda_matrix(self, lambd):
        size = len(self.matrix)
        return Matrix([[(lambd if i == j else 0) - self.matrix[i][j] for j in range(size)] for i in range(size)])

    # Đưa về ma trận tối giản dòng
    def row_reduced_echelon_form(self):
        A = [row[:] for row in self.matrix]  # Sao chép ma trận
        n = len(A)
        m = len(A[0])

        lead = 0
        for r in range(n):
            if lead >= m:
                return A
            i = r
            while A[i][lead] == 0:
                i += 1
                if i == n:
                    i = r
                    lead += 1
                    if m == lead:
                        return A
            A[i], A[r] = A[r], A[i]
            lv = A[r][lead]
            A[r] = [element / lv for element in A[r]]
            for i in range(r + 1, n):
                lv = A[i][lead]
                A[i] = [A[i][j] - lv * A[r][j] for j in range(m)]
            lead += 1
        return A

    def find_eigenvector(self, eigenvalue):
        lambda_matrix = self.create_lambda_matrix(eigenvalue)

        rref_matrix = lambda_matrix.row_reduced_echelon_form() # ma trận đã ở dạng bậc thang

        # Tạo vector riêng từ ma trận RREF

        # vector 0 luôn là nghiệm tầm thường
        n = len(self.matrix)
        eigenvector = [0] * n

        # Tìm vector riêng từ ma trận RREF
        # Bước này sẽ trả về một vector không tầm thường nếu ma trận có nghiệm

        # xác định các dòng 0 trong ma trận bậc thang
        free_vars = [i for i in range(n) if all(rref_matrix[j][i] == 0 for j in range(n))]

        # ở mỗi dòng 0 của ma trận (lambda*I_n - A_n) đã về dạng bậc thang, v_i eigenvector tương ứng có thể bằng bất kì giá trị nào, ta lấy giá trị 1 cho v_i luôn
        for i in free_vars:
            eigenvector[i] = 1

        # những vị trí v còn lại của tọa độ eigenvector này vẫn bằng 0
        # nếu ko có dòng 0 nào thì eigenvector duy nhất vẫn là vector 0 (nghiệm tầm thường)

        return eigenvector

    def eigenvectors(self, eigenvalues):
        vectors = []
        for e in eigenvalues:
            vectors.append(self.find_eigenvector(e))
        return vectors

    # Cach 2: cho ma tran kich thuoc lon
    # Power Iteration: Tim gia tri rieng lon nhat
    # Later...

    # - Chéo hóa ma trận
    def diagonalize(self, eigenvalues, eigenvectors):

        if len(eigenvalues) < self.C:
            raise ValueError("Không thể chéo hóa ma trận vì số lượng giá trị riêng nhỏ hơn kích thước ma trận!")

        if len(eigenvectors) < self.C:
            raise ValueError("Không thể chéo hóa ma trận vì số lượng vector riêng nhỏ hơn kích thước ma trận!")

        # Tạo ma trận P từ các eigenvector
        P = Matrix(eigenvectors).transpose()

        # Tạo ma trận đường chéo D từ các eigenvalue
        D = Matrix([[0] * self.C for _ in range(self.C)])

        for i, eigenvalue in enumerate(eigenvalues):
            D.matrix[i][i] = eigenvalue

        return P, D

    # - Kiếm tra xem, ma trận có xác định dương hay không
    def check_positive_define(self, eigenvalues):
        return all(v > 0 for v in eigenvalues)

    # - Tính chuẩn của 1 vector và ma trận.

    # Chuan l2 cua vector
    @staticmethod
    def vector_norm_l2(vector):
        return math.sqrt(sum(x ** 2 for x in vector))

    # Chuan l1 cua vector
    @staticmethod
    def vector_norm_l1(vector):
        return sum(abs(x) for x in vector)

    # Chuan fr0benius cua ma tran
    def matrix_norm_frobenius(self):
        return math.sqrt(sum(cell ** 2 for row in self.matrix for cell in row))

    def shape(self):
        return self.R, self.C  # so hang, so cot

    def __str__(self):
        return '\n'.join([' '.join(map(str, row)) for row in self.matrix])
