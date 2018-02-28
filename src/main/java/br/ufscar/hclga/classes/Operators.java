package br.ufscar.hclga.classes;

import java.util.Random;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Operators {

//Use these if initiate the rules using seeding
//To include relation tests, put one more "<=" in the end of numericOperators
private static String[] numericOperatorsRel = {"<=", ">=", "<= <=", "<="};
private static String[] numericOperatorsNRel = {"<=", ">=", "<= <="};

private static String[] categoricOperators = {"=", "!=", "in"};

/* ===========================================================
* Randomly select a numeric operator
* =========================================================== */
public static double getInitialFlagValue(double probabilityFlag) {
    double flagValue = 0.0;

    Random generator = new Random();
    double num = generator.nextDouble();

    if (num <= probabilityFlag) {
        flagValue = 1.0;
    }

    return flagValue;
}

/* ===========================================================
* Randomly select a numeric operator
* =========================================================== */
public static double getNumericOperator() {
    double operator;

    Random generator = new Random();
    operator = generator.nextInt(numericOperatorsNRel.length);

    return operator;
}

/* ===========================================================
* Randomly select a categoric operator
* =========================================================== */
public static double getCategoricOperator() {
    Random generator = new Random();
    double operator = generator.nextInt(categoricOperators.length);

    return operator;
}

}
