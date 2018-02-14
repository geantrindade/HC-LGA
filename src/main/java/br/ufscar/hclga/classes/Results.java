package br.ufscar.hclga.classes;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Results {

    /* ===========================================================
     * Mount a string with the condition of a given term
     * =========================================================== */
    public static ArrayList<String> getStringConditionValues(double[] rule, int indexTerm) {

        ArrayList<String> termValues = new ArrayList<String>();

        int posAttribute = indexTerm / 4;
        double numericOperatorValue = rule[indexTerm + 1];

        if (Datasets.getInfoAttributes().get(posAttribute) == 1) {
            //Numeric attribute
            //0 -> "<="
            //1 -> ">="
            //2 -> "<= <="
            //3 -> attribute <= attribute

            switch ((int) numericOperatorValue) {
                case 0:
                    termValues.add(" <= " + rule[indexTerm + 3]);
                    break;

                case 1:
                    termValues.add(" >= " + rule[indexTerm + 2]);
                    break;

                case 2:
                    termValues.add(rule[indexTerm + 2] + " <= ");
                    termValues.add(" <= " + rule[indexTerm + 3]);
                    break;

                case 3:
                    termValues.add(" <= ");
                    break;
            }
        } else {
            //Categoric attribute
            //0 -> "="
            //1 -> "!="
            //2 -> "in"

            int posCatAttribute = Datasets.getCategoricAttributes().indexOf(posAttribute);
            double infLim = 0;
            String categoricValue = "";

            switch ((int) numericOperatorValue) {
                case 0:
                    infLim = rule[indexTerm + 2];
                    categoricValue = Datasets.getCategoricMapping().get(posCatAttribute).get((int) infLim)[0];
                    termValues.add(" = " + categoricValue);
                    break;

                case 1:
                    infLim = rule[indexTerm + 2];
                    categoricValue = Datasets.getCategoricMapping().get(posCatAttribute).get((int) infLim)[0];
                    termValues.add(" != " + categoricValue);
                    break;

                case 2:
                    infLim = rule[indexTerm + 2];
                    String[] values = Datasets.getCategoricMapping().get(posCatAttribute).get((int) infLim);
                    categoricValue = mountCompoundStringCatValues(values);
                    termValues.add(" in " + categoricValue);
                    break;
            }
        }

        return termValues;
    }

    /* ===========================================================
     * Mount a string with possible values for the "in" operator
     * =========================================================== */
    private static String mountCompoundStringCatValues(String[] values) {

        String stringValues = "";

        stringValues += "{";

        for (int i = 0; i < (values.length - 1); i++) {
            stringValues += values[i] + ", ";
        }

        stringValues += values[values.length - 1] + "}";

        return stringValues;
    }

    /* ===========================================================
     * Returns the complete term of the rule
     * =========================================================== */
    private static String mountCompleteTerm(double[] rule, int indexTerm, ArrayList<String> termValues) {

        String ruleTerm = "";

        int posAttribute = indexTerm / 4;
        double numericOperatorValue = rule[indexTerm + 1];
        ArrayList<String> namesAttributes = Datasets.getNamesAttributes();

        if (Datasets.getInfoAttributes().get(posAttribute) == 1 && numericOperatorValue == 2) {

            ruleTerm = termValues.get(0) + namesAttributes.get(posAttribute) + termValues.get(1);

        } else if (Datasets.getInfoAttributes().get(posAttribute) == 1 && numericOperatorValue == 3) {

            ruleTerm = namesAttributes.get((int) rule[indexTerm + 2]) + termValues.get(0) + namesAttributes.get((int) rule[indexTerm + 3]);

        } else {
            ruleTerm = namesAttributes.get(posAttribute) + termValues.get(0);
        }

        return ruleTerm;
    }

    /* ===========================================================
     * Print the rules
     * =========================================================== */
    public static void printRules(ArrayList<Individual> population) {

        for (int i = 0; i < population.size(); i++) {
            double[] rule = population.get(i).getRule();
            ArrayList<Integer> activeTerms = GeneticOperators.getPosActiveTerms(rule);
            System.out.print("Rule " + i + ": ");
            for (int j = 0; j < activeTerms.size(); j++) {
                int indexTerm = activeTerms.get(j);

                //Returns a string with the condition of a given term
                ArrayList<String> termValues = Results.getStringConditionValues(rule, indexTerm);

                //Returns the complete term of the rule
                String ruleTerm = mountCompleteTerm(rule, indexTerm, termValues);

                System.out.print(ruleTerm);
                if (j < activeTerms.size() - 1) {
                    System.out.print(" AND ");
                }
            }
            System.out.println();
        }
    }

    public static void printRule(double[] rule) {

        ArrayList<Integer> activeTerms = GeneticOperators.getPosActiveTerms(rule);

        System.out.print("Rule = ");
        for (int j = 0; j < activeTerms.size(); j++) {
            int indexTerm = activeTerms.get(j);

            //Returns a string with the condition of a given term
            ArrayList<String> termValues = Results.getStringConditionValues(rule, indexTerm);

            //Returns the complete term of the rule
            String ruleTerm = mountCompleteTerm(rule, indexTerm, termValues);

            System.out.print(ruleTerm);
            if (j < activeTerms.size() - 1) {
                System.out.print(" AND ");
            }
        }
        System.out.println();
    }

    /* ===========================================================
     * Gets a string corresponding to a rule
     * =========================================================== */
    public static String getRule(double[] rule) {

        ArrayList<Integer> activeTerms = GeneticOperators.getPosActiveTerms(rule);
        String stringRule = "";

        int indexTerm;
        ArrayList<String> termValues;
        String ruleTerm;

        for (int j = 0; j < activeTerms.size(); j++) {
            indexTerm = activeTerms.get(j);

            //Returns a string with the condition of a given term
            termValues = Results.getStringConditionValues(rule, indexTerm);

            //Returns the complete term of the rule
            ruleTerm = mountCompleteTerm(rule, indexTerm, termValues);

            stringRule += ruleTerm;

            if (j < activeTerms.size() - 1) {
                stringRule += " AND ";
            }
        }

        return stringRule;

    }

    /* ===========================================================
     * Obtain the consequent of a rule
     * =========================================================== */
 /*public static Individual obtainConsequent(Individual individual) {
    
     int numCoveredExamples = individual.getNumberCoveredExamples();
     ArrayList<Integer> indexCoveredExamples = individual.getIndexCoveredExamples();
     double[] consequentVector = individual.getConsequentVector();
    
     for (int i = 0; i < indexCoveredExamples.size(); i++) {
    
     int indexExample = indexCoveredExamples.get(i);
     int[] binaryClasses = Classes.getBinaryClassesTrain().get(indexExample);
    
     for (int j = 0; j < binaryClasses.length; j++) {
     consequentVector[j] += binaryClasses[j];
     }
     }
    
     for (int i = 0; i < consequentVector.length; i++) {
     consequentVector[i] = consequentVector[i] / numCoveredExamples;
     }
    
     individual.setConsequentVector(consequentVector);
    
     return individual;
     }
     */
 /* ===========================================================
     * Save the prediction obtained
     * =========================================================== */
    public static void savePredictions(double[][] matrixPredictions, double AUPRCvalue, ArrayList<double[]> AUPRCClasses) {

        try {

            String filePredictions = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/" + Paths.getNameFilePredictions();

            String fileResultsAUPRCClasses = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/" + Paths.getNameFileAUPRCclasses();

            //Prediction real numbers
            File predictions = new File(filePredictions);
            FileWriter fstream = new FileWriter(predictions);
            BufferedWriter out = new BufferedWriter(fstream);

            for (int i = 0; i < matrixPredictions.length; i++) {
                for (int j = 0; j < matrixPredictions[i].length; j++) {
                    out.write(matrixPredictions[i][j] + " ");
                }
                out.write("\n");
            }
            //AUPRC overall
            out.write("\nAU(PRC) value = " + AUPRCvalue);

            out.close();

            //Results in all classes
            File resultsClasses = new File(fileResultsAUPRCClasses);
            FileWriter fstreamRClasses = new FileWriter(resultsClasses);
            BufferedWriter outRClasses = new BufferedWriter(fstreamRClasses);

            for (int ind = 0; ind < Classes.getClasses().length; ind++) {

                outRClasses.write(Classes.getClasses()[ind] + ", ");
                outRClasses.write("AUPRC = " + AUPRCClasses.get(ind)[0] + "\n");
            }

            outRClasses.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /* ===========================================================
     * Save the fmeasure per level (Single-Label version)
     * =========================================================== */
    public static void saveFmeasureRun(double[] fmeasureLevels) {

        try {

            String fileFmeasureLevels = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/" + Paths.getNameFileFmeasureLevels();

            //Results in each level
            File fmeasures = new File(fileFmeasureLevels);
            FileWriter fstream = new FileWriter(fmeasures);
            BufferedWriter out = new BufferedWriter(fstream);
            int numLevel = 0;

            for (int level = 0; level < Parameters.getNumLevels(); level++) {
                numLevel++;
                out.write("Fmeasure level " + numLevel + " = " + fmeasureLevels[level] + "\n");
            }

            out.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /* ===========================================================
     * Save the training times of a run
     * =========================================================== */
    public static void saveTrainingTimesRun(double[] trainingTimes) {

        try {

            String fileTrainingTimesRun = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/" + Paths.getNameFileTrainingTimes();

            //Results in each level
            File tTimes = new File(fileTrainingTimesRun);
            FileWriter fstream = new FileWriter(tTimes);
            BufferedWriter out = new BufferedWriter(fstream);

            out.write("--------------- Training time ---------------\n");
            out.write("Milliseconds = " + trainingTimes[0] + "\n");
            out.write("Seconds = " + trainingTimes[1] + "\n");
            out.write("Minutes = " + trainingTimes[2] + "\n");
            out.write("Hours = " + trainingTimes[3] + "\n");

            out.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /* ===========================================================
     * Save the mean training times of all runs
     * =========================================================== */
    public static void saveMeanTrainingTimes(double[] meanMsTimes, double[] meanSTimes,
            double[] meanMTimes, double[] meanHTimes) {

        try {

            String fileMeanTrainingTimes = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/" + Paths.getNameFileMeanTrainingTimes();

            //Results in each level
            File tTimes = new File(fileMeanTrainingTimes);
            FileWriter fstream = new FileWriter(tTimes);
            BufferedWriter out = new BufferedWriter(fstream);

            out.write("--------------- Training time ---------------\n");
            out.write("Milliseconds = " + meanMsTimes[0] + "(" + meanMsTimes[1] + ")\n");
            out.write("Seconds = " + meanSTimes[0] + "(" + meanSTimes[1] + ")\n");
            out.write("Minutes = " + meanMTimes[0] + "(" + meanMTimes[1] + ")\n");
            out.write("Hours = " + meanHTimes[0] + "(" + meanHTimes[1] + ")\n");

            out.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /* ===========================================================
     * Save the prediction obtained in a given threshold
     * =========================================================== */
    public static void savePredictionsThreshold(int[][] matrixPredictions, double[] evalResults,
            ArrayList<double[]> evalResultsClasses, double threshold) {

        try {

            String filePredictions = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/Thresholds/" + threshold + "/" + Paths.getNameFilePredictions();

            String filePRclasses = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/Thresholds/" + threshold + "/" + Paths.getNameFilePRclasses();

            //Predictions
            File predictions = new File(filePredictions);
            FileWriter fstream = new FileWriter(predictions);
            BufferedWriter out = new BufferedWriter(fstream);

            for (int i = 0; i < matrixPredictions.length; i++) {
                String predClasses = "";
                if (Parameters.getHierarchyType().equals("Tree")) {
                    predClasses = retrievePredictionTree(matrixPredictions[i], Classes.getClasses());
                } else {
                    //predClasses = retrievePredictionDAG(matrixPredictions[i], Classes.getClasses());
                }

                out.write(predClasses + "\n");
            }
            //Precision and Recall values
            out.write("\n\n");
            out.write("Precision = " + evalResults[0] + '\n');
            out.write("Recall = " + evalResults[1] + '\n');

            out.write("\tTrue Positives = " + evalResults[2] + '\n');
            out.write("\tFalse Positives = " + evalResults[3] + '\n');
            out.write("\tTotal Real = " + evalResults[4] + '\n');

            out.close();

            //Precision and recall for each class
            File PRclasses = new File(filePRclasses);
            FileWriter fstreamPRclasses = new FileWriter(PRclasses);
            BufferedWriter outPRclasses = new BufferedWriter(fstreamPRclasses);

            double[] results;

            for (int i = 0; i < Classes.getClasses().length; i++) {

                results = evalResultsClasses.get(i);

                outPRclasses.write(Classes.getClasses()[i] + ": ");

                outPRclasses.write("Precision = " + results[0] + ", ");
                outPRclasses.write("Recall = " + results[1] + ", ");

                outPRclasses.write("\tTrue Positives = " + results[2] + ", ");
                outPRclasses.write("\tFalse Positives = " + results[3] + ", ");
                outPRclasses.write("\tTotal Real = " + results[4] + "\n");

            }

            outPRclasses.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /*===========================================================================
     * Retrieve the predictions for a tree structure
     *===========================================================================*/
    static String retrievePredictionTree(int[] vectorOutputs,
            String[] classes) {

        ArrayList<String> predictedClasses = new ArrayList<String>();

        for (int index = (vectorOutputs.length - 1); index >= 0; index--) {
            if (vectorOutputs[index] == 1) {
                if (index == (vectorOutputs.length - 1)) {
                    predictedClasses.add(classes[index]);
                } else {
                    int index2 = predictedClasses.size();
                    String nameClass = classes[index];
                    nameClass = "^" + nameClass + "/";
                    Pattern pattern = Pattern.compile(nameClass);
                    int found = 0;
                    Matcher m;

                    for (int k = 0; k < index2; k++) {
                        m = pattern.matcher(predictedClasses.get(k));
                        if (m.find()) {
                            found = 1;
                            break;
                        }
                    }
                    if (found == 0) {
                        predictedClasses.add(0, classes[index]);
                    }
                }
            }
        }

        String predClasses = "";
        for (int i = 0; i < predictedClasses.size(); i++) {
            if (i == 0) {
                predClasses = predictedClasses.get(i);
            } else {
                predClasses = predClasses + "@" + predictedClasses.get(i);
            }
        }
        return predClasses;
    }

    /* ===========================================================
     * Save the rules obtained
     * =========================================================== */
    public static void saveRules(ArrayList<Individual> individuals) {

        try {

            String fileRules = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/" + Paths.getNameFileRules();

            //Prediction real numbers
            File rules = new File(fileRules);
            FileWriter fstream = new FileWriter(rules);
            BufferedWriter out = new BufferedWriter(fstream);

            double[] rule;
            String stringRule;
            double[] consequent;

            for (int i = 0; i < individuals.size(); i++) {

                rule = individuals.get(i).getRule();
                stringRule = getRule(rule);

                out.write("Rule " + i + " = " + stringRule + " THEN ");

                consequent = individuals.get(i).getMeanClassLabelVectorCovered();

                if (Parameters.getMultiLabel() == 0) {
                    consequent = getHigherProbabilities(consequent);
                }

                for (int j = 0; j < consequent.length; j++) {
                    out.write(consequent[j] + " ");
                }
                out.write("\n\n");
            }
            out.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

      public static void saveRules(ArrayList<Individual> individuals, String finalPath) {

        try {

            //Prediction real numbers
            File rules = new File(finalPath + Paths.getNameFileRules());
            FileWriter fstream = new FileWriter(rules);
            BufferedWriter out = new BufferedWriter(fstream);

            double[] rule;
            String stringRule;
            double[] consequent;

            for (int i = 0; i < individuals.size(); i++) {

                rule = individuals.get(i).getRule();
                stringRule = getRule(rule);

                out.write("Rule " + i + " = " + stringRule + " THEN ");

                consequent = individuals.get(i).getMeanClassLabelVectorCovered();

                if (Parameters.getMultiLabel() == 0) {
                    consequent = getHigherProbabilities(consequent);
                }

                for (int j = 0; j < consequent.length; j++) {
                    out.write(consequent[j] + " ");
                }
                out.write("\n\n");
            }
            out.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

  
    
    /* ===========================================================
     * Obtain the predictions
     * =========================================================== */
    public static double[][] obtainPredictions(ArrayList<Individual> individuals, double[] defaultRule) {

        double[][] matrixPredictions = new double[Classes.getBinaryClassesTest().size()][Classes.getBinaryClassesTest().get(0).length];

        ArrayList<String[]> datasetTest = Datasets.getDatasetTest();

        String[] example;
        int coverage;

        for (int i = 0; i < datasetTest.size(); i++) {

            example = datasetTest.get(i);
            coverage = 0;

            Individual individual;

            //Verify if the example is covered by a rule
            for (int j = 0; j < individuals.size(); j++) {
                individual = individuals.get(j);
                coverage = individual.verifyRule(example, individual.getRule(), 0);

                if (coverage == 1) {
                    double[] vectorConsequent = individual.getMeanClassLabelVectorCovered();
                    System.arraycopy(vectorConsequent, 0, matrixPredictions[i], 0, vectorConsequent.length);
                    break;
                }
            }

            //If no rule classify the example, apply the default rule
            if (coverage == 0) {
                System.arraycopy(defaultRule, 0, matrixPredictions[i], 0, defaultRule.length);
            }
        }

        return matrixPredictions;
    }

    /*===========================================================================
     * Save the interpolated points for individual classes (For testing)
     *===========================================================================*/
    static void saveInterpolationClasses(ArrayList<ArrayList<ArrayList<Double>>> interpolatedPrecisionPoints,
            ArrayList<ArrayList<ArrayList<Double>>> interpolatedRecallPoints) {

        String fileInterpolationsClasses = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/" + Paths.getNameFileInterpolationClasses();

        try {

            File interpolated = new File(fileInterpolationsClasses);
            FileWriter FI = new FileWriter(interpolated);
            BufferedWriter outFI = new BufferedWriter(FI);

            //Iterates over all classes
            for (int ind = 0; ind < Classes.getClasses().length; ind++) {

                ArrayList<ArrayList<Double>> precision = interpolatedPrecisionPoints.get(ind);
                ArrayList<ArrayList<Double>> recall = interpolatedRecallPoints.get(ind);

                outFI.write("================== ");
                outFI.write(Classes.getClasses()[ind]);
                outFI.write(" ==================\n");

                outFI.write("Precision\t\t\t\t\tRecall\n");

                for (int i = 0; i < recall.size(); i++) {
                    for (int j = 0; j < recall.get(i).size(); j++) {

                        outFI.write(precision.get(i).get(j).toString());
                        outFI.write("\t\t\t\t\t" + recall.get(i).get(j).toString() + "\n");
                    }
                }
            }

            outFI.write("================== END ==================\n");

            outFI.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /*===========================================================================
     * Save the interpolated points (For testing)
     *===========================================================================*/
    static void saveInterpolation(ArrayList<ArrayList<Double>> recall,
            ArrayList<ArrayList<Double>> precision) {

        String fileInterpolations = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                + "/" + Parameters.getFileDatasetTest() + "/Run" + HC_LGA.getNumRun() + "/" + Paths.getNameFileInterpolation();

        try {

            File interpolated = new File(fileInterpolations);
            FileWriter FI = new FileWriter(interpolated);
            BufferedWriter outFI = new BufferedWriter(FI);

            outFI.write("Precision\t\t\t\t\tRecall\n");

            for (int i = 0; i < recall.size(); i++) {
                for (int j = 0; j < recall.get(i).size(); j++) {

                    outFI.write(precision.get(i).get(j).toString());
                    outFI.write("\t\t\t\t\t" + recall.get(i).get(j).toString() + "\n");
                }
            }

            outFI.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /*===========================================================================
     * Calculate mean and standard deviation of a vector of results
     *===========================================================================*/
    static double[] calculateMeanSd(double[] vector) {

        double[] meanSd = new double[3];
        double sum = 0;
        double mean = 0;
        double sd = 0;

        for (int i = 0; i < vector.length; i++) {
            sum += vector[i];
        }
        mean = sum / (vector.length);
        meanSd[0] = mean;
        sum = 0;

        for (int i = 0; i < vector.length; i++) {
            sum += Math.pow((vector[i] - meanSd[0]), 2);
        }
        sum = sum / (vector.length - 1);
        sd = Math.sqrt(sum);

        meanSd[1] = sd;

        return meanSd;
    }

    /*===========================================================================
     * Calculate mean and standard deviation for the f-measures in each level
     *===========================================================================*/
    static ArrayList<double[]> calculateMeanSdFmesureLevels(ArrayList<double[]> fmeasures) {

        ArrayList<double[]> meanSd = new ArrayList<double[]>();

        for (int level = 0; level < Parameters.getNumLevels(); level++) {

            double[] meanSdLevel = new double[2];
            double sum = 0;
            double mean = 0;
            double sd = 0;

            for (int i = 0; i < fmeasures.size(); i++) {
                sum += fmeasures.get(i)[level];
            }

            mean = sum / (fmeasures.size());
            meanSdLevel[0] = mean;
            sum = 0;

            for (int i = 0; i < fmeasures.size(); i++) {
                sum += Math.pow((fmeasures.get(i)[level] - meanSdLevel[0]), 2);
            }
            sum = sum / (fmeasures.size() - 1);
            sd = Math.sqrt(sum);

            meanSdLevel[1] = sd;

            meanSd.add(meanSdLevel);
        }

        return meanSd;
    }

    /*===========================================================================
     * Calculate the means of the AUPRC of all executions for all classes 
     *===========================================================================*/
    public static ArrayList<double[]> calculateMeansAUPRCClasses(ArrayList<ArrayList<double[]>> allAUPRCClasses) {

        ArrayList<double[]> meansAllAUPRCClasses = new ArrayList<double[]>();

        //Store the values of all executions for one class
        double[] AUPRCClass = new double[allAUPRCClasses.size()];

        //Go through all classes
        for (int i = 0; i < allAUPRCClasses.get(0).size(); i++) {

            //Go through all executions
            for (int j = 0; j < allAUPRCClasses.size(); j++) {

                AUPRCClass[j] = allAUPRCClasses.get(j).get(i)[0];

            }

            //Calculate the mean and sd
            double[] meanSd = calculateMeanSd(AUPRCClass);
            meanSd[2] = allAUPRCClasses.get(0).get(i)[1];
            meansAllAUPRCClasses.add(meanSd);

        }

        return meansAllAUPRCClasses;
    }

    /*===========================================================================
     * Save the means of AUPRC of all classes
     *===========================================================================*/
    static void saveMeansAUPRCClasses(ArrayList<double[]> meansAllAUPRCClasses,
            double[] meanAUPRC, double[] allAUPRCTest, int sizeDataset) {

        String fileMeanAUPRCclasses = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                + "/" + Parameters.getFileDatasetTest() + "/" + Paths.getNameFileMeanAUPRCclasses();

        try {
            //Results in all classes and overall
            File resultsClasses = new File(fileMeanAUPRCclasses);
            FileWriter fstreamRClasses = new FileWriter(resultsClasses);
            BufferedWriter outRClasses = new BufferedWriter(fstreamRClasses);
            int numRun = 0;

            for (int i = 0; i < allAUPRCTest.length; i++) {
                numRun++;
                outRClasses.write("AUPRC Run " + numRun + " = " + allAUPRCTest[i] + "\n");
            }

            outRClasses.write("===================================================================\n");
            outRClasses.write("Mean (Sd) test AU.PRC = " + meanAUPRC[0] + " (" + meanAUPRC[1] + ")\n");
            outRClasses.write("===================================================================\n\n");

            /*outR.write("----------------------- Mean training time -----------------------\n");
             outR.write("Mean (Sd) milliseconds = " + meanMsTimes[0] + " (" + meanMsTimes[1] + ")\n");
             outR.write("Mean (Sd) seconds = " + meanSTimes[0] + " (" + meanSTimes[1] + ")\n");
             outR.write("Mean (Sd) minutes = " + meanMTimes[0] + " (" + meanMTimes[1] + ")\n");
             outR.write("Mean (Sd) hours = " + meanHTimes[0] + " (" + meanHTimes[1] + ")\n");*/
            for (int ind = 0; ind < Classes.getClasses().length; ind++) {

                double freqClass = meansAllAUPRCClasses.get(ind)[2] / sizeDataset;

                outRClasses.write(Classes.getClasses()[ind] + ", ");
                outRClasses.write("Mean AUPRC = " + meansAllAUPRCClasses.get(ind)[0]);
                outRClasses.write(" (" + meansAllAUPRCClasses.get(ind)[1] + ")");
                outRClasses.write(", Freq:" + freqClass + "\n");

            }

            outRClasses.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

    /*===========================================================================
     * For single-label, lets get only the higher probabilities
     * in the vectors for each level
     *===========================================================================*/
    public static double[][] getHigherProbabilities(double[][] matrixPredictions) {

        String[] classes = Classes.getClasses();

        double[][] matrixPredictionsSingleLabel = new double[matrixPredictions.length][classes.length];
        /*for(int i=0; i<matrixPredictionsSingleLabel.length; i++){
         for(int j=0; j<matrixPredictionsSingleLabel[i].length; j++){
         matrixPredictionsSingleLabel[i][j] = -1;
         }
         }*/

        String patternClass;
        Pattern pattern;
        double higher;
        int indexHigher;
        Matcher m;

        for (int i = 0; i < matrixPredictions.length; i++) {
            patternClass = "";

            for (int level = 0; level < Parameters.getNumLevels(); level++) {
                pattern = Pattern.compile("^" + patternClass + "[0-9]+$");

                higher = 0;
                indexHigher = 0;

                for (int indexClass = 0; indexClass < classes.length; indexClass++) {
                    m = pattern.matcher(classes[indexClass]);

                    if (m.find()) {
                        if (matrixPredictions[i][indexClass] > higher) {
                            higher = matrixPredictions[i][indexClass];
                            indexHigher = indexClass;
                        }
                    }
                }

                if (higher > 0) {
                    matrixPredictionsSingleLabel[i][indexHigher] = higher;
                    patternClass = classes[indexHigher] + "/";
                } else {
                    break;
                }
            }
        }

        return matrixPredictionsSingleLabel;
    }

    public static double[] getHigherProbabilities(double[] vectorPredictions) {

        String[] classes = Classes.getClasses();

        double[] vectorPredictionsSingleLabel = new double[classes.length];

        String patternClass = "";
        Pattern pattern;
        double higher;
        int indexHigher;
        Matcher m;

        for (int level = 0; level < Parameters.getNumLevels(); level++) {
            pattern = Pattern.compile("^" + patternClass + "[0-9]+$");

            higher = 0;
            indexHigher = 0;

            for (int indexClass = 0; indexClass < classes.length; indexClass++) {
                m = pattern.matcher(classes[indexClass]);

                if (m.find()) {
                    if (vectorPredictions[indexClass] > higher) {
                        higher = vectorPredictions[indexClass];
                        indexHigher = indexClass;
                    }
                }
            }

            if (higher > 0) {
                vectorPredictionsSingleLabel[indexHigher] = higher;
                patternClass = classes[indexHigher] + "/";
            } else {
                break;
            }
        }

        return vectorPredictionsSingleLabel;
    }

    /*===========================================================================
     * Save the mean f-measure values for each level (Single-Label version)
     *===========================================================================*/
    public static void saveMeanFmeasureLevels(ArrayList<double[]> meanFmeasures) {

        try {

            String fileMeanFmeasureLevels = Paths.getPathResults() + "/" + Parameters.getHierarchyType()
                    + "/" + Parameters.getFileDatasetTest() + "/" + Paths.getNameFileMeanFmeasureLevels();

            //Results in each level
            File fmeasures = new File(fileMeanFmeasureLevels);
            FileWriter fstream = new FileWriter(fmeasures);
            BufferedWriter out = new BufferedWriter(fstream);
            int numLevel = 0;

            for (int level = 0; level < Parameters.getNumLevels(); level++) {
                numLevel++;
                out.write("Mean Fmeasure level " + numLevel + " = " + meanFmeasures.get(level)[0]);
                out.write(" (" + meanFmeasures.get(level)[1] + ")\n");
            }

            out.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }
}
