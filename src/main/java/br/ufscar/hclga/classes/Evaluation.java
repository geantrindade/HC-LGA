package br.ufscar.hclga.classes;

import java.util.ArrayList;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Evaluation {

private double AUPRC;
private ArrayList<double[]> AUPRCClasses;
private double[] fmeasureLevels;

/*===========================================================================
* Evaluation method to calculate the AU(PRC) values
*===========================================================================*/
public void evaluationAUPRC(double[][] matrixPredictions) {
    //Store precision and recall values
    ArrayList<double[]> valuesPrecisionRecall = new ArrayList<double[]>();
    ArrayList<ArrayList<double[]>> valuesPrecisionRecallClasses = new ArrayList<ArrayList<double[]>>();
    ArrayList<Double> thresholdValues = Parameters.getThresholdValues();

    //Iterate over all thresholds
    for (int indexThres = 0; indexThres < thresholdValues.size(); indexThres++) {
        //Matrix to store the outputs on the test data after applying thresholds
        int[][] binaryMatrix = new int[Classes.getBinaryClassesTest().size()][Classes.getBinaryClassesTest().get(0).length];

        //Threshold values used
        double threshold = thresholdValues.get(indexThres) / 100;

        //Apply the threshold
        applyThresholds(binaryMatrix, matrixPredictions, threshold, 1);

        //Hierarchical Precision and Recall evaluation metrics
        double[] evalResults = evaluationPrecRec(Classes.getBinaryClassesTest(), binaryMatrix);
        valuesPrecisionRecall.add(evalResults);

        //Hierarchical Precision and Recall evaluation for each class individually
        ArrayList<double[]> evalResultsClasses = evaluationPrecRecClasses(Classes.getBinaryClassesTest(), binaryMatrix);
        valuesPrecisionRecallClasses.add(evalResultsClasses);

        //Save the predictions and results for each threshold
        //Results.savePredictionsThreshold(binaryMatrix, evalResults, evalResultsClasses, thresholdValues.get(indexThres));
    }
    //Calculate AU(PRC)
    AUPRC = calculateAUPRC(valuesPrecisionRecall);

    //Calculate AU(PRC) of each class
    AUPRCClasses = calculateAUPRCClasses(valuesPrecisionRecallClasses);
}

/*===========================================================================
* Evaluation method to calculate the AU(PRC) values
*===========================================================================*/
public double evaluationAUPRCFitness(double[][] matrixPredictions, ArrayList<Integer> indexExamples) {
    //Store precision and recall values
    ArrayList<double[]> valuesPrecisionRecall = new ArrayList<double[]>();
    ArrayList<Double> thresholdValues = Parameters.getThresholdValues();

    double AUPRCFitness = 0;

    //Iterate over all thresholds
    for (int indexThres = 0; indexThres < thresholdValues.size(); indexThres++) {
        //Matrix to store the outputs on the test data after applying thresholds
        int[][] binaryMatrix = new int[indexExamples.size()][Classes.getClasses().length];

        //Threshold values used
        double threshold = thresholdValues.get(indexThres) / 100;

        //Apply the threshold
        applyThresholds(binaryMatrix, matrixPredictions, threshold, 0);

        ArrayList<int[]> trueClasses = new ArrayList<int[]>();
        int[] classes = new int[Classes.getClasses().length];

        for (int i = 0; i < indexExamples.size(); i++) {
            int posExample = indexExamples.get(i);
            System.arraycopy(Classes.getBinaryClassesTrain().get(posExample), 0, classes, 0, classes.length);
            trueClasses.add(classes);
        }

        //Hierarchical Precision and Recall evaluation metrics
        double[] evalResults = evaluationPrecRec(trueClasses, binaryMatrix);
        valuesPrecisionRecall.add(evalResults);

    }
    //Calculate AU(PRC)
    AUPRCFitness = calculateAUPRCFitness(valuesPrecisionRecall);

    return AUPRCFitness;
}

/*===========================================================================
* Evaluation measure to calculated the Fmeasure per level (Single-Label)
* Will consider that if the value assign to a class is higher then 0, the
* class is predicted (1). Will not penalize over-specialized predictions
*===========================================================================*/
public void evaluationFmeasure(double[][] matrixPredictions) {
    fmeasureLevels = new double[Parameters.getNumLevels()];

    //Matrix to store the outputs on the test data
    int[][] binaryMatrix = new int[Classes.getBinaryClassesTest().size()][Classes.getBinaryClassesTest().get(0).length];

    applyThresholds(binaryMatrix, matrixPredictions, 0, 0);

    //Get position of classes by level
    ArrayList<ArrayList<Integer>> positionClassesLevel = Classes.getPositionClassesLevel();

    //Store the positions of the classes that will be used in the evaluation
    //ArrayList<Integer> positionsClassesToEvaluate = new ArrayList<Integer>();
    for (int i = 0; i < Parameters.getNumLevels(); i++) {
        //positionsClassesToEvaluate.addAll(positionClassesLevel.get(i));
        //F-measure per level
        double fmeasure = evaluationFmeasureLevel(Classes.getBinaryClassesTest(), binaryMatrix, positionClassesLevel.get(i));
        fmeasureLevels[i] = fmeasure;
    }
}

/*===========================================================================
* Calculate the F-measure for a level, given an array list with the 
* position of the classes for the level
* Will not penalize over-specialization
*===========================================================================*/
static double evaluationFmeasureLevel(ArrayList<int[]> trueClasses, int[][] predictedClasses, ArrayList<Integer> positionClassesLevel) {
    double fmeasure = 0;
    double sumIntersection = 0;
    double minSumPredicted = 0;
    double sumReal = 0;

    for (int numInst = 0; numInst < trueClasses.size(); numInst++) {
        double sumPredictedExample = 0;
        double sumRealExample = 0;

        for (int i = 0; i < positionClassesLevel.size(); i++) {
            int posClass = positionClassesLevel.get(i);

            if (predictedClasses[numInst][posClass] == 1 && trueClasses.get(numInst)[posClass] == 1) {
                sumIntersection++;
            }

            if (predictedClasses[numInst][posClass] == 1) {
                sumPredictedExample++;
            }

            if (trueClasses.get(numInst)[posClass] == 1) {
                sumRealExample++;
                sumReal++;
            }
        }

        //Get the minimum value. This will not penalize over-specialization
        if (sumPredictedExample < sumRealExample) {
            minSumPredicted += sumPredictedExample;

        } else {
            minSumPredicted += sumRealExample;
        }
    }

    //Hierarchical Precision
    double hPrecision = 0.0;

    if (minSumPredicted != 0) {
        hPrecision = sumIntersection / minSumPredicted;
    }

    //Hierarchical Recall
    double hRecall = 0.0;

    if (sumReal != 0) {
        hRecall = sumIntersection / sumReal;
    }

    //Fmeasure
    if (hPrecision != 0 || hRecall != 0) {
        fmeasure = (2 * hPrecision * hRecall) / (hPrecision + hRecall);
    }

    return fmeasure;
}

/*===========================================================================
* Calculate the F-measure for a level, given an array list with the 
* position of the classes for the level
* Will not penalize over-specialization
*===========================================================================*/
public double evaluationFmeasureFitness(double[][] predictedClasses, ArrayList<Integer> indexExamples) {
    double fmeasure = 0;
    double sumIntersection = 0;
    double minSumPredicted = 0;
    double sumReal = 0;

    //Matrix to store the outputs on the test data
    int[][] binaryMatrix = new int[indexExamples.size()][Classes.getClasses().length];

    applyThresholds(binaryMatrix, predictedClasses, 0, 0);

    for (int i = 0; i < indexExamples.size(); i++) {
        int numInst = indexExamples.get(i);
        double sumPredictedExample = 0;
        double sumRealExample = 0;

        for (int j = 0; j < Classes.getClasses().length; j++) {
            if (binaryMatrix[i][j] == 1 && Classes.getBinaryClassesTrain().get(numInst)[j] == 1) {
                sumIntersection++;
            }

            if (binaryMatrix[i][j] == 1) {
                sumPredictedExample++;
            }

            if (Classes.getBinaryClassesTrain().get(numInst)[j] == 1) {
                sumRealExample++;
                sumReal++;
            }
        }

        //Get the minimum value. This will not penalize over-specialization
        if (sumPredictedExample < sumRealExample) {
            minSumPredicted += sumPredictedExample;

        } else {
            minSumPredicted += sumRealExample;
        }
    }

    //Hierarchical Precision
    double hPrecision = 0.0;

    if (minSumPredicted != 0) {
        hPrecision = sumIntersection / minSumPredicted;
    }

    //Hierarchical Recall
    double hRecall = 0.0;

    if (sumReal != 0) {
        hRecall = sumIntersection / sumReal;
    }

    //Fmeasure
    if (hPrecision != 0 || hRecall != 0) {
        fmeasure = (2 * hPrecision * hRecall) / (hPrecision + hRecall);
    }

    return fmeasure;
}

public double evaluationRecallFitness(double[][] predictedClasses, ArrayList<Integer> indexExamples) {
    double sumIntersection = 0;
    double sumReal = 0;

    //Matrix to store the outputs on the test data
    int[][] binaryMatrix = new int[indexExamples.size()][Classes.getClasses().length];

    applyThresholds(binaryMatrix, predictedClasses, 0, 0);

    for (int i = 0; i < indexExamples.size(); i++) {
        int numInst = indexExamples.get(i);

        for (int j = 0; j < Classes.getClasses().length; j++) {
            if (binaryMatrix[i][j] == 1 && Classes.getBinaryClassesTrain().get(numInst)[j] == 1) {
                sumIntersection++;
            }

            if (Classes.getBinaryClassesTrain().get(numInst)[j] == 1) {
                sumReal++;
            }
        }
    }

    //Hierarchical Recall
    double hRecall = sumIntersection / sumReal;

    return hRecall;
}

public double evaluationPrecisionFitness(double[][] predictedClasses, ArrayList<Integer> indexExamples) {
    double sumIntersection = 0;
    double sumPredicted = 0;

    //Matrix to store the outputs on the test data
    int[][] binaryMatrix = new int[indexExamples.size()][Classes.getClasses().length];

    applyThresholds(binaryMatrix, predictedClasses, 0, 0);

    for (int i = 0; i < indexExamples.size(); i++) {
        int numInst = indexExamples.get(i);

        for (int j = 0; j < Classes.getClasses().length; j++) {
            if (binaryMatrix[i][j] == 1 && Classes.getBinaryClassesTrain().get(numInst)[j] == 1) {
                sumIntersection++;
            }

            if (binaryMatrix[i][j] == 1) {
                sumPredicted++;
            }
        }
    }

    //Hierarchical Precision
    double hPrecision = sumIntersection / sumPredicted;

    return hPrecision;
}

/*===========================================================================
* Calculate the F-measure for a given prediction of an example
* Will not penalize over-specialization
*===========================================================================*/
public double fMeasurePrediction(int indexExample, double[] prediction) {
    double fmeasure = 0;
    double sumIntersection = 0;
    double minSumPredicted = 0;
    double sumReal = 0;
    double sumPredicted = 0;

    for (int j = 0; j < Classes.getClasses().length; j++) {
        if (prediction[j] > 0 && Classes.getBinaryClassesTrain().get(indexExample)[j] == 1) {
            sumIntersection++;
        }

        if (prediction[j] > 0) {
            sumPredicted++;
        }

        if (Classes.getBinaryClassesTrain().get(indexExample)[j] == 1) {
            sumReal++;
        }
    }

    //Get the minimum value. This will not penalize over-specialization
    if (sumPredicted < sumReal) {
        minSumPredicted = sumPredicted;

    } else {
        minSumPredicted = sumReal;
    }

    //Hierarchical Precision
    double hPrecision = 0.0;

    if (minSumPredicted != 0) {
        hPrecision = sumIntersection / minSumPredicted;
    }

    //Hierarchical Recall
    double hRecall = 0.0;

    if (sumReal != 0) {
        hRecall = sumIntersection / sumReal;
    }

    //Fmeasure
    if (hPrecision != 0 || hRecall != 0) {
        fmeasure = (2 * hPrecision * hRecall) / (hPrecision + hRecall);
    }

    return fmeasure;
}

/*===========================================================================
* Calculate the AU(PRC) for classes individually
*===========================================================================*/
static ArrayList<double[]> calculateAUPRCClasses(ArrayList<ArrayList<double[]>> valuesPrecisionRecallClasses) {
    ArrayList<double[]> AUPRCClasses = new ArrayList<double[]>();
    ArrayList<ArrayList<ArrayList<Double>>> interpolatedPrecisionPoints = new ArrayList<ArrayList<ArrayList<Double>>>();
    ArrayList<ArrayList<ArrayList<Double>>> interpolatedRecallPoints = new ArrayList<ArrayList<ArrayList<Double>>>();

    //Iterates over all classes
    int aux = 0;

    for (int ind = 0; ind < Classes.getClasses().length; ind++) {
        ArrayList<ArrayList<Double>> precision = new ArrayList<ArrayList<Double>>();
        ArrayList<ArrayList<Double>> recall = new ArrayList<ArrayList<Double>>();
        double[] AUPRC = new double[2];
        int count = 0;
        double totalTrueClass = 0;

        for (int i = valuesPrecisionRecallClasses.size() - 1; i > 0; i--) {
            //Recover data for interpolation
            double[] dataInterpolation = getDataInterpolation(valuesPrecisionRecallClasses.get(i).get(aux),
                    valuesPrecisionRecallClasses.get(i - 1).get(aux));

            totalTrueClass = dataInterpolation[9];

            //Get points between A and B to interpolate
            ArrayList<ArrayList<Double>> points = getPoints(dataInterpolation, count);

            if (i < (valuesPrecisionRecallClasses.size() - 1)) {
                precision.get(precision.size() - 1).remove(precision.get(precision.size() - 1).size() - 1);
                recall.get(recall.size() - 1).remove(recall.get(recall.size() - 1).size() - 1);
            }

            precision.add(points.get(0));
            recall.add(points.get(1));
            count++;
        }

        interpolatedPrecisionPoints.add(precision);
        interpolatedRecallPoints.add(recall);

        AUPRC[0] = calculateAreaUnderCurve(recall, precision);
        AUPRC[1] = totalTrueClass;

        AUPRCClasses.add(AUPRC);
        aux++;
    }

    Results.saveInterpolationClasses(interpolatedPrecisionPoints, interpolatedRecallPoints);

    return AUPRCClasses;
}

/*===========================================================================
* Calculate the AU(PRC)
*===========================================================================*/
static double calculateAUPRC(ArrayList<double[]> valuesPrecisionRecall) {
    double AUPRC = 0;
    ArrayList<ArrayList<Double>> precision = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Double>> recall = new ArrayList<ArrayList<Double>>();
    int count = 0;

    for (int i = valuesPrecisionRecall.size() - 1; i > 0; i--) {
        //Recover data for interpolation
        double[] dataInterpolation = getDataInterpolation(valuesPrecisionRecall.get(i), valuesPrecisionRecall.get(i - 1));

        //Get points between A and B to interpolate
        ArrayList<ArrayList<Double>> points = getPoints(dataInterpolation, count);

        if (i < (valuesPrecisionRecall.size() - 1)) {
            precision.get(precision.size() - 1).remove(precision.get(precision.size() - 1).size() - 1);
            recall.get(recall.size() - 1).remove(recall.get(recall.size() - 1).size() - 1);
        }

        precision.add(points.get(0));
        recall.add(points.get(1));
        count++;
    }

    AUPRC = calculateAreaUnderCurve(recall, precision);

    Results.saveInterpolation(recall, precision);

    return AUPRC;
}

static double calculateAUPRCFitness(ArrayList<double[]> valuesPrecisionRecall) {
    double AUPRC = 0;
    ArrayList<ArrayList<Double>> precision = new ArrayList<ArrayList<Double>>();
    ArrayList<ArrayList<Double>> recall = new ArrayList<ArrayList<Double>>();
    int count = 0;

    for (int i = valuesPrecisionRecall.size() - 1; i > 0; i--) {
        //Recover data for interpolation
        double[] dataInterpolation = getDataInterpolation(valuesPrecisionRecall.get(i), valuesPrecisionRecall.get(i - 1));

        //Get points between A and B to interpolate
        ArrayList<ArrayList<Double>> points = getPoints(dataInterpolation, count);

        if (i < (valuesPrecisionRecall.size() - 1)) {
            precision.get(precision.size() - 1).remove(precision.get(precision.size() - 1).size() - 1);
            recall.get(recall.size() - 1).remove(recall.get(recall.size() - 1).size() - 1);
        }

        precision.add(points.get(0));
        recall.add(points.get(1));
        count++;
    }

    AUPRC = calculateAreaUnderCurve(recall, precision);

    return AUPRC;
}

/*===========================================================================
* Calculate the area under a curve
*===========================================================================*/
static double calculateAreaUnderCurve(ArrayList<ArrayList<Double>> recall,
        ArrayList<ArrayList<Double>> precision) {
    double AUPRC = 0;
    ArrayList<Double> x = new ArrayList<Double>();
    ArrayList<Double> y = new ArrayList<Double>();

    for (int i = 0; i < recall.size(); i++) {
        for (int j = 0; j < recall.get(i).size(); j++) {
            x.add(recall.get(i).get(j));
            y.add(precision.get(i).get(j));
        }
    }

    for (int i = 0; i < x.size() - 1; i++) {
        AUPRC += (x.get(i + 1) - x.get(i)) * y.get(i + 1)
                + (x.get(i + 1) - x.get(i)) * (y.get(i) - y.get(i + 1)) / 2;
    }

    return AUPRC;
}

/*===========================================================================
* Interpolate points between two P/R values
*===========================================================================*/
static ArrayList<ArrayList<Double>> getPoints(double[] dataInterpolation, int count) {
    ArrayList<ArrayList<Double>> points = new ArrayList<ArrayList<Double>>();
    ArrayList<Double> prec = new ArrayList<Double>();
    ArrayList<Double> reca = new ArrayList<Double>();

    double localSkew = dataInterpolation[0];

    double tpA = dataInterpolation[1];
    double fpA = dataInterpolation[2];
    double tpB = dataInterpolation[3];
    double fpB = dataInterpolation[4];

    double precA = dataInterpolation[5];
    double recaA = dataInterpolation[6];
    double precB = dataInterpolation[7];
    double recaB = dataInterpolation[8];

    double total = dataInterpolation[9];

    double param = tpB - tpA;

    double newPrec;
    double newReca;

    if (count == 0 && recaA > 0) {
        prec.add(precA);
        reca.add(0.0);
    }

    if (count == 0) {
        prec.add(precA);
        reca.add(recaA);
    }

    if (count > 0 && (precA != 0 || recaA != 0 || precB != 0 || recaB != 0)) {
        prec.add(precA);
        reca.add(recaA);
    }

    for (int tp = (int) tpA + 1; tp < tpB; tp++) {
        double fp = fpA + localSkew * (tp - tpA);
        newPrec = tp / (tp + fp);
        newReca = tp / total;

        prec.add(newPrec);
        reca.add(newReca);
    }

    prec.add(precB);
    reca.add(recaB);

    points.add(prec);
    points.add(reca);

    return points;
}

/*===========================================================================
* Get data values to interpolate two PR points
*===========================================================================*/
static double[] getDataInterpolation(double[] valuesPRA, double[] valuesPRB) {
    double[] dataInterpolation = new double[10];
    double localSkew;

    //Values for point A
    double precisionA = valuesPRA[0];
    double recallA = valuesPRA[1];
    double tpA = valuesPRA[2];
    double fpA = valuesPRA[3];
    double total = valuesPRA[4];

    //Values for point B
    double precisionB = valuesPRB[0];
    double recallB = valuesPRB[1];
    double tpB = valuesPRB[2];
    double fpB = valuesPRB[3];

    if ((tpB - tpA) == 0) {
        localSkew = 0;

    } else {
        localSkew = (fpB - fpA) / (tpB - tpA);
    }

    dataInterpolation[0] = localSkew;
    dataInterpolation[1] = tpA;
    dataInterpolation[2] = fpA;
    dataInterpolation[3] = tpB;
    dataInterpolation[4] = fpB;
    dataInterpolation[5] = precisionA;
    dataInterpolation[6] = recallA;
    dataInterpolation[7] = precisionB;
    dataInterpolation[8] = recallB;
    dataInterpolation[9] = total;

    return dataInterpolation;
}

/* ===========================================================
* Apply a threshold to the results obtained
* =========================================================== */
static void applyThresholds(int[][] binaryMatrix, double[][] matrixPredictions, double threshold, int multilabel) {
    for (int numInst = 0; numInst < binaryMatrix.length; numInst++) {
        for (int i = 0; i < binaryMatrix[0].length; i++) {
            binaryMatrix[numInst][i] = getOutputThreshold(matrixPredictions[numInst][i], threshold, multilabel);
        }
    }
}

/*===========================================================================
* Gets a binary prediction given a vector of real values and a threshold
*===========================================================================*/
static int getOutputThreshold(double realValue, double threshold, int multilabel) {
    int output = 0;

    if (multilabel == 1) {
        if (realValue >= threshold) {
            output = 1;

        } else {
            output = 0;
        }

    } else if (realValue > threshold) {
        output = 1;

    } else {
        output = 0;
    }

    return output;
}

/*===========================================================================
* Hierarchical Precision and Recall evaluation metrics
*===========================================================================*/
static double[] evaluationPrecRec(ArrayList<int[]> trueClasses, int[][] predictedClasses) {
    //Store the results
    double[] evalResults = new double[5];

    //Sum of predicted and real classes
    double sumIntersection = 0;
    double sumPredicted = 0;
    double sumReal = 0;
    double FP = 0;

    if (Parameters.getHierarchyType().equals("Tree")) {
        for (int numInst = 0; numInst < trueClasses.size(); numInst++) {
            for (int i = 0; i < trueClasses.get(0).length; i++) {
                if (predictedClasses[numInst][i] == 1 && trueClasses.get(numInst)[i] == 1) {
                    sumIntersection++;
                }

                if (predictedClasses[numInst][i] == 1) {
                    sumPredicted++;
                }

                if (predictedClasses[numInst][i] == 1 && trueClasses.get(numInst)[i] == 0) {
                    FP++;
                }

                if (trueClasses.get(numInst)[i] == 1) {
                    sumReal++;
                }
            }
        }

    } else {
        String[] illegalGOclasses = Classes.getIllegalGOclasses();

        //Get positions of the illegal classes
        ArrayList<Integer> illegalPositions = new ArrayList<Integer>();

        for (int pos = 0; pos < illegalGOclasses.length; pos++) {
            for (int pos2 = 0; pos2 < Classes.getClasses().length; pos2++) {
                if (illegalGOclasses[pos].equals(Classes.getClasses()[pos2])) {
                    illegalPositions.add(pos2);
                    break;
                }
            }
        }

        for (int numInst = 0; numInst < trueClasses.size(); numInst++) {
            for (int i = 0; i < trueClasses.get(0).length; i++) {
                if (illegalPositions.contains(i) == false) {
                    if (predictedClasses[numInst][i] == 1 && trueClasses.get(numInst)[i] == 1) {
                        sumIntersection++;
                    }

                    if (predictedClasses[numInst][i] == 1) {
                        sumPredicted++;
                    }

                    if (predictedClasses[numInst][i] == 1 && trueClasses.get(numInst)[i] == 0) {
                        FP++;
                    }

                    if (trueClasses.get(numInst)[i] == 1) {
                        sumReal++;
                    }
                }
            }
        }
    }

    //Hierarchical Precision
    double hPrecision = 0.0;

    if (sumPredicted != 0) {
        hPrecision = sumIntersection / sumPredicted;
    }

    //Hierarchical Recall
    double hRecall = 0.0;

    if (sumReal != 0) {
        hRecall = sumIntersection / sumReal;
    }

    evalResults[0] = hPrecision;
    evalResults[1] = hRecall;
    evalResults[2] = sumIntersection; //TP
    evalResults[3] = FP;              //FP
    evalResults[4] = sumReal;         //True

    return evalResults;
}

/*===========================================================================
* Hierarchical Precision and Recall evaluation metrics for each class
*===========================================================================*/
static ArrayList<double[]> evaluationPrecRecClasses(ArrayList<int[]> trueClasses, int[][] predictedClasses) {
    //Store the results
    ArrayList<double[]> evalResultsClasses = new ArrayList<double[]>();

    //Iterates over all classes
    for (int i = 0; i < trueClasses.get(0).length; i++) {
        double[] evalResults = new double[5];

        //Sum of predicted and real classes
        double sumIntersection = 0;
        double sumPredicted = 0;
        double sumReal = 0;
        double FP = 0;

        //Iterates over all test instances
        for (int numInst = 0; numInst < trueClasses.size(); numInst++) {
            if (predictedClasses[numInst][i] == 1 && trueClasses.get(numInst)[i] == 1) {
                sumIntersection++;
            }

            if (predictedClasses[numInst][i] == 1) {
                sumPredicted++;
            }

            if (predictedClasses[numInst][i] == 1 && trueClasses.get(numInst)[i] == 0) {
                FP++;
            }

            if (trueClasses.get(numInst)[i] == 1) {
                sumReal++;
            }
        }

        //Hierarchical Precision
        double hPrecision = 0.0;

        if (sumPredicted != 0) {
            hPrecision = sumIntersection / sumPredicted;
        }

        //Hierarchical Recall
        double hRecall = 0.0;

        if (sumReal != 0) {
            hRecall = sumIntersection / sumReal;
        }

        evalResults[0] = hPrecision;
        evalResults[1] = hRecall;
        evalResults[2] = sumIntersection; //TP
        evalResults[3] = FP;              //FP
        evalResults[4] = sumReal;         //True

        evalResultsClasses.add(evalResults);
    }

    return evalResultsClasses;
}

public double getAUPRC() {
    return AUPRC;
}

public double[] getFmeasureLevels() {
    return fmeasureLevels;
}

public ArrayList<double[]> getAUPRCClasses() {
    return AUPRCClasses;
}

}
