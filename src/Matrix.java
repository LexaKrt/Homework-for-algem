import java.util.Scanner;

public class Matrix {
    private final int rowsCount;
    private final int columnCount;

    private final double[][] matrix;

    private int numForDeterminant = 1;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
        this.rowsCount = matrix.length;
        this.columnCount = matrix[0].length;
    }


    public Matrix kramerMethod(Matrix freeCoefs) {
        Matrix result = new Matrix(new double[freeCoefs.rowsCount][freeCoefs.columnCount]);
        if (freeCoefs.rowsCount != columnCount | freeCoefs.columnCount != 1) {
            System.out.println("Error! Wrong data! Zero matrix has been returned");
            return result;
        }
        double mainDeterminant = determinant();
        for (int i = 0; i < columnCount; i++) {   // columnCount == rowsCount
            Matrix tempMatrix = copyMatrix();
            tempMatrix.swapColumnsForCramer(i, freeCoefs);
            result.matrix[i][0] = round(tempMatrix.determinant() / mainDeterminant);
        }
        return result;
    }


    private Matrix copyMatrix() {
        Matrix result = new Matrix(new double[rowsCount][columnCount]);
        for (int i = 0; i < rowsCount; i++) {
            for (int j = 0; j < columnCount; j++) {
                result.matrix[i][j] = matrix[i][j];
            }
        }
        return result;
    }


    private void swapColumnsForCramer(int indOfFirstMatrix, Matrix secondMatrix) {
        for (int i = 0; i < rowsCount; i++) {
            matrix[i][indOfFirstMatrix] = secondMatrix.matrix[i][0];
        }
    }

    public Matrix solve(Matrix freeCoefficients) {
        if (freeCoefficients.matrix[0].length != 1 | freeCoefficients.matrix.length != rowsCount) {
            System.out.println("Wrong matrix-column!");
            return null;
        }
        if (determinant() == 0) {
            System.out.println("Matrix can't be reversed!");
            return null;
        } else {
            Matrix reversedMatrix = createOnesMatrix();
            gaussMethodForSolve(reversedMatrix);
            reversedGaussMethodForSolve(reversedMatrix);
            for (int i = 0; i < rowsCount; i++) {
                double coefficient = matrix[i][i];
                for (int j = 0; j < columnCount; j++) {
                    reversedMatrix.matrix[i][j] /= coefficient;
                    reversedMatrix.matrix[i][j] = round(reversedMatrix.matrix[i][j]);
                }
            }
            return reversedMatrix.multiplicationOfMatrices(freeCoefficients);
        }
    }


    public void gaussMethodForSolve(Matrix anotherMatrix) {
        for (int i = 0; i < rowsCount - 1; i++) {
            if (matrix[i][i] == 0) {
                int indexThatNotZero = indexNotZeroInColumn(i);
                if (indexThatNotZero != i) {
                    swapLines(i, indexThatNotZero);
                    anotherMatrix.swapLines(i, indexThatNotZero);
                } else {
                    indexThatNotZero = indexNotZeroInRow(i);
                    if (indexThatNotZero != i) {
                        swapColumns(i, indexThatNotZero);
                        anotherMatrix.swapColumns(i, indexThatNotZero);
                    } else {
                        i++;
                        continue;
                    }
                }
            }
            for (int m = i + 1; m < rowsCount; m++) {
                double k = - matrix[m][i] / matrix[i][i];
                for (int j = 0; j < columnCount; j++) {
                    matrix[m][j] += matrix[i][j] * k;
                    anotherMatrix.matrix[m][j] += anotherMatrix.matrix[i][j] * k;
                    matrix[m][j] = round(matrix[m][j]);
                    anotherMatrix.matrix[m][j] = round(anotherMatrix.matrix[m][j]);
                }
            }
        }
    }

    public void reversedGaussMethodForSolve(Matrix anotherMatrix) {
        for (int i = rowsCount - 1; i > 0; i--) {
            for (int m = i - 1; m > -1; m--) {
                double k = - matrix[m][i] / matrix[i][i];
                for (int j = 0; j < columnCount; j++) {
                    matrix[m][j] += matrix[i][j] * k;
                    anotherMatrix.matrix[m][j] += anotherMatrix.matrix[i][j] * k;
                    matrix[m][j] = round(matrix[m][j]);
                    anotherMatrix.matrix[m][j] = round(anotherMatrix.matrix[m][j]);
                }
            }
        }
    }


    private Matrix createOnesMatrix() {
        Matrix result = new Matrix(new double[rowsCount][columnCount]);
        for (int i = 0; i < rowsCount; i++) {
            result.matrix[i][i] = 1;
        }
        return result;
    }


    private boolean isAllZeros() {
        for (int i = 0; i < rowsCount; i++) {
            for (int j = 0; j < columnCount; j++) {
                if (matrix[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    public double determinant() {
        Matrix result = copyMatrix();
        if (result.rowsCount != result.columnCount) {
            return 0;
        }
        result.gaussMethod();
        double determinant = 1;
        for (int i = 0; i < rowsCount; i++) {
            determinant *= result.matrix[i][i];
        }
        return round(determinant * result.numForDeterminant);
    }


    public void reversedGaussMethod() {
        for (int i = rowsCount - 1; i > 0; i--) {
            for (int m = i - 1; m > -1; m--) {
                double k = - matrix[m][i] / matrix[i][i];
                for (int j = 0; j < columnCount; j++) {
                    matrix[m][j] += matrix[i][j] * k;
                    matrix[m][j] = round(matrix[m][j]);
                }
            }
        }
    }


    public double round(double num) {
        return (double) Math.round(num * 1000000) / 1000000;
    }

//    rank
    public int rank() {
        gaussMethod();
        int rank = rowsCount;

        for (int i = 0; i < rowsCount; i++) {
            boolean flag = true;
            for (int j = 0; j < columnCount; j++) {
                if (matrix[i][j] != 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                rank--;
            }
        }
        return rank;
    }

// block of gaussMethod
    public void gaussMethod() {
        for (int i = 0; i < rowsCount - 1; i++) {
            if (matrix[i][i] == 0) {
                int indexThatNotZero = indexNotZeroInColumn(i);
                if (indexThatNotZero != i) {
                    swapLines(i, indexThatNotZero);
                } else {
                    indexThatNotZero = indexNotZeroInRow(i);
                    if (indexThatNotZero != i) {
                        swapColumns(i, indexThatNotZero);
                    } else {
                        i++;
                        continue;
                    }
                }
            }
             for (int m = i + 1; m < rowsCount; m++) {
                double k = - matrix[m][i] / matrix[i][i];
                for (int j = i; j < columnCount; j++) {
                    matrix[m][j] += matrix[i][j] * k;
                    matrix[m][j] = round(matrix[m][j]);
                }
            }
        }
    }


    private int indexNotZeroInRow(int indWhereZero) {
        for (int i = indWhereZero + 1; i < columnCount; i++) {
            if (matrix[indWhereZero][i] != 0) {
                return i;
            }
        }
        return indWhereZero;
    }


    private void swapColumns(int indWhereZero, int indNotZero) {
        for (int i = 0; i < rowsCount; i++) {
            double temp = matrix[i][indWhereZero];
            matrix[i][indWhereZero] = matrix[i][indNotZero];
            matrix[i][indNotZero] = temp;
        }
        numForDeterminant *= -1;
    }


    private int indexNotZeroInColumn(int indWhereZero) {
        for (int i = indWhereZero + 1; i < rowsCount; i++) {
            if (matrix[i][indWhereZero] != 0) {
                return i;
            }
        }
        return indWhereZero;
    }


    private void swapLines(int indWhereZero, int indNotZero) {
            double[] temp = matrix[indWhereZero];
            matrix[indWhereZero] = matrix[indNotZero];
            matrix[indNotZero] = temp;
            numForDeterminant *= -1;
    }


    private void swapLines(int indWhereZero) {
        for (int i = indWhereZero; i < rowsCount; i++) {
            if (matrix[i][indWhereZero] != 0) {
                double[] temp = matrix[indWhereZero];
                matrix[indWhereZero] = matrix[i];
                matrix[i] = temp;
                numForDeterminant *= -1;
                break;
            }
        }
    }
//      end of gaussMethod block

    public Matrix multiplicationOfMatrices(Matrix anotherMatrix) {
        Scanner scanner = new Scanner(System.in);
        if (this.columnCount == anotherMatrix.rowsCount) {
            Matrix result = new Matrix(new double[this.rowsCount][anotherMatrix.columnCount]);
            cycleForMultiplication(result, anotherMatrix);
            return result;
        } else if (anotherMatrix.columnCount == this.rowsCount) {
            System.out.println("We can't multiply matrix A on matrix B.");
            System.out.print("Do you want to multiply second matrix on first matrix?: ");
            String answer = scanner.next();
            if (answer.equalsIgnoreCase("yes") || answer.equalsIgnoreCase("да") || answer.equalsIgnoreCase("д") || answer.equalsIgnoreCase("y")) {
                Matrix result = new Matrix(new double[anotherMatrix.rowsCount][this.columnCount]);
                anotherMatrix.cycleForMultiplication(result, this);
                return result;
            } else {
                System.out.println("Null has been returned.");
                return null;
            }
        } else {
            System.out.println("We can't multiply matrices!");
            return null;
        }
    }


    private void cycleForMultiplication(Matrix result, Matrix anotherMatrix) {
        for (int i = 0; i < this.rowsCount; i++) {
            for (int j = 0; j < anotherMatrix.columnCount; j++) {
                for (int k = 0; k < anotherMatrix.rowsCount; k++) {
                    result.matrix[i][j] += matrix[i][k] * anotherMatrix.matrix[k][j];
                    result.matrix[i][j] = round(result.matrix[i][j]);
                }
            }
        }
    }


    public Matrix sumOfMatrices(Matrix anotherMatrix) {
        if ((this.rowsCount) == (anotherMatrix.rowsCount) & (this.columnCount) == (anotherMatrix.columnCount)) {
            Matrix result = new Matrix(new double[rowsCount][columnCount]);
            for (int i = 0; i < this.rowsCount; i++) {
                for (int j = 0; j < this.columnCount; j++) {
                     result.matrix[i][j] = this.matrix[i][j] + anotherMatrix.matrix[i][j];
                }
            }
            return result;
        } else {
            System.out.println("We can't summary these Matrices!");
            return null;
        }
    }


    public void fillMatrixElemByElem() {
        Scanner scanner = new Scanner(System.in);
        printMatrix();
        for (int i = 0; i < rowsCount; i++) {
            for (int j = 0; j < columnCount; j++) {
                System.out.printf("Enter the value of %d row and %d column: ", i + 1, j + 1);
                String lineWithValue = scanner.nextLine();
                System.out.println();
                matrix[i][j] = Double.parseDouble(lineWithValue);
                printMatrix();
            }
        }
    }


    public void fillMatrixRowByRow() {
        Scanner scanner = new Scanner(System.in);
        printMatrix();
        for (int i = 0; i < rowsCount; i++) {
            System.out.printf("Enter the %d row with \" \"(space) symbol: ", i + 1);
            String lineWithValue = scanner.nextLine();
            String[] arr = lineWithValue.split(" ");
            System.out.println();

            while (arr.length != columnCount) {
                System.out.println("Wrong data! Try again.");
                System.out.printf("Enter the %d row with \" \"(space) symbol: ", i + 1);
                lineWithValue = scanner.nextLine();
                arr = lineWithValue.split(" ");

            }

            for (int j = 0; j < columnCount; j++) {
                matrix[i][j] = Double.parseDouble(arr[j]);
            }

            printMatrix();
        }
    }


    public void printMatrix() {
        System.out.printf("Your matrix %dx%d: ", rowsCount, columnCount);
        System.out.println();
        for (int i = 0;i < rowsCount; i++) {
            for (int j = 0;j < columnCount; j++) {
                System.out.print(matrix[i][j] + "  ");
            }
            System.out.println();
        }
        System.out.println();
    }


    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        String RED = "\u001B[38;2;255;0;0m";
        String RESET = "\u001B[0m";


        String answer = "passWord";
        while (!answer.equalsIgnoreCase("start")) {
            System.out.print("Enter \"start\" to run a program: ");
            answer = scanner.next();
            System.out.println();
        }

        System.out.println(RED + "\tMatrix methods  " + RESET);

        /*Matrix matrix1 = new Matrix(setSizeOfNewMatrix());
        matrix1.chooseMethodToFillMatrix();

        Matrix matrix2 = new Matrix(setSizeOfNewMatrix());
        matrix2.chooseMethodToFillMatrix();

        Matrix resultOfMultiplication = matrix1.multiplicationOfMatrices(matrix2);*/



       /* Matrix resultOfSummarize = matrix1.sumOfMatrices(matrix2);

        System.out.println("This is result of summarize of two matrices: ");
        resultOfSummarize.printMatrix();*/



        /*System.out.println("This is result of multiplying of two matrices");
        resultOfMultiplication.printMatrix();*/


        /*Matrix matrix3 = new Matrix(setSizeOfNewMatrix());
        matrix3.chooseMethodToFillMatrix();


        System.out.println("Rank of your matrix is: " + matrix3.rank());*/



       /* Matrix matrix4 = new Matrix(matrix3.matrix);
        matrix4.gaussMethod();

        matrix4.reversedGaussMethod();

        matrix4.printMatrix();*/

        double[][] array = {
                {0, 2, 0},
                {0, 0, 1},
                {1, 0, 0}
        };

        double[][] freeCoefs = {
                {1}, {2}, {3}
        };

        /*Matrix matrix5 = new Matrix(setSizeOfNewMatrix());
        matrix5.chooseMethodToFillMatrix();

        System.out.println(RED + "Enter matrix of free coefficients" + RESET);
        Matrix freeCoefficients = new Matrix(setSizeOfNewMatrix());
        freeCoefficients.chooseMethodToFillMatrix();

        Matrix resultSolve = matrix5.solve(freeCoefficients);
        resultSolve.printMatrix();*/



        /*Matrix matrix6 = new Matrix(setSizeOfNewMatrix());
        matrix6.chooseMethodToFillMatrix();

        System.out.println(RED + "Enter matrix of free coefficients" + RESET);
        Matrix freeCoefs = new Matrix(setSizeOfNewMatrix());
        freeCoefs.chooseMethodToFillMatrix();

        Matrix resultKramer = matrix6.kramerMethod(freeCoefs);
        resultKramer.printMatrix();*/
    }


    private static double[][] setSizeOfNewMatrix() {
        Scanner scanner = new Scanner(System.in);
        int tempRowsNumber;
        int tempColumnsNumber;

        System.out.println("Enter the number of rows and number of columns: ");
        System.out.print("\t-number of rows: ");
        tempRowsNumber = scanner.nextInt();
        System.out.println();

        while (tempRowsNumber < 1) {
            System.out.println("Wrong data! Try again.");
            System.out.println();
            System.out.println("\t-number of rows: ");
            tempRowsNumber = scanner.nextInt();
        }

        System.out.print("\t-number of Columns: ");
        tempColumnsNumber = scanner.nextInt();
        System.out.println();

        while (tempColumnsNumber < 1) {
            System.out.println("Wrong data! Try again.");
            System.out.println();
            System.out.println("\t-number of Columns: ");
            tempColumnsNumber = scanner.nextInt();
        }

        return new double[tempRowsNumber][tempColumnsNumber];
    }


    public void chooseMethodToFillMatrix() {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Enter the number of method to fill matrix (1 or 2): ");
        System.out.println("1) Fill matrix element by element");
        System.out.print("2) Fill matrix row by row" + "\n" + "Your answer is: ");
        int answer = scanner.nextInt();
        System.out.println();

        while (answer != 1 & answer != 2) {
            System.out.println("Wrong data! Try again.");
            System.out.println();
            System.out.println("Enter the number of method to fill matrix (1 or 2): ");
            System.out.println("1) Fill matrix element by element");
            System.out.print("2) Fill matrix row by row" + "\n" + "Your answer is: ");
            answer = scanner.nextInt();
        }

        if (answer == 1) {
            fillMatrixElemByElem();
        } else {
            fillMatrixRowByRow();
        }
    }
}