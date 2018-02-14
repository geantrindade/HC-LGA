package br.ufscar.hclga.classes;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Datasets {

private static ArrayList<String[]> datasetTrain;
private static ArrayList<String[]> datasetValid;
private static ArrayList<String[]> datasetTest;
private static ArrayList<Integer> infoAttributes;
private static ArrayList<String> namesAttributes;
private static ArrayList<Integer> categoricAttributes;
private static ArrayList<String> categoricValues;
private static ArrayList<ArrayList<String[]>> categoricMapping;
private ArrayList<String[]> possibleCombinations;
private String tokenAttribute = "@ATTRIBUTE";
private String tokenData = "@DATA";
private String tokenNumeric = "numeric";
private String tokenHierarchical = "hierarchical";
private static String tokenMissingValue = "?";
private static String tokenRootClass = "root";
private static ArrayList<double[]> meanValues = new ArrayList<double[]>();
private static ArrayList<ArrayList<String>> modeValues = new ArrayList<ArrayList<String>>();

public Datasets() {
    datasetTrain = new ArrayList<String[]>();
//        datasetValid = new ArrayList<String[]>();

    String pathDatasetTrain = Paths.getPathDatasets() + Parameters.getFileDatasetTrain();
//        String pathDatasetValid = Paths.getPathDatasets() + Parameters.getFileDatasetValid();

    //Will store info about the attribute types, numeric or categoric
    infoAttributes = new ArrayList<Integer>();

    //Names of the attributes
    namesAttributes = new ArrayList<String>();

    //Information about the categorical attributes
    categoricAttributes = new ArrayList<Integer>();
    categoricValues = new ArrayList<String>();
    categoricMapping = new ArrayList<ArrayList<String[]>>();

    //Number of the attribute
    int numAttribute = -1;

    //Compile the tokes to search in the files
    Pattern patternAttribute = Pattern.compile(tokenAttribute, Pattern.CASE_INSENSITIVE);
    Pattern patternData = Pattern.compile(tokenData, Pattern.CASE_INSENSITIVE);
    Pattern patternNumeric = Pattern.compile(tokenNumeric, Pattern.CASE_INSENSITIVE);
    Pattern patternHierarchical = Pattern.compile(tokenHierarchical, Pattern.CASE_INSENSITIVE);

    try {
        //Read training dataset and get infos about attribute values
        //--> 1 == Numeric Attribute / 2 == Categoric Attribute
        //==================================================================
        FileReader readerTrain = new FileReader(pathDatasetTrain);
        BufferedReader rTrain = new BufferedReader(readerTrain);

        String line = null;
        int dataFound = 0;

        while ((line = rTrain.readLine()) != null) {
            //Get attributes information util do not find @DATA token
            if (dataFound == 0) {
                Matcher mAttribute = patternAttribute.matcher(line);
                Matcher mData = patternData.matcher(line);

                //See if reached @DATA token
                if (mData.find()) {
                    dataFound = 1;
                } //See if @ATTRIBUTE token was found
                else if (mAttribute.find()) {
                    numAttribute++;

                    //If so, check if attribute is hierarchical, numeric or categoric
                    Matcher mNumeric = patternNumeric.matcher(line);
                    Matcher mHierarchical = patternHierarchical.matcher(line);

                    if (mHierarchical.find()) {
                        //Build structure to store the classes
                        //Tree or DAG
                        if ("Tree".equals(Parameters.getHierarchyType())) {
                            Classes.setTreeClasses(line, tokenHierarchical);

                        } else if ("DAG".equals(Parameters.getHierarchyType())) {
                            Classes.setDAGclasses(line, tokenHierarchical);
                        }

                    } else if (mNumeric.find()) {
                        infoAttributes.add(1);

                        //Save the attribute name
                        String attrName = getAttributeName(line, 1);
                        namesAttributes.add(attrName);

                    } else {
                        infoAttributes.add(2);

                        //Save the attribute name
                        String attrName = getAttributeName(line, 2);
                        namesAttributes.add(attrName);

                        //Store the possible values for the categoric attribute
                        String stringValues = getCategoricValues(line);
                        categoricAttributes.add(numAttribute);
                        categoricValues.add(stringValues);
                        ArrayList<String[]> mapping = mapPossibleCategoricValues(stringValues);
                        categoricMapping.add(mapping);
                    }
                }
            } //Token @DATA was found
            else {
                String[] vetLine = line.split(",");
                datasetTrain.add(vetLine);
            }
        }

        rTrain.close();
        readerTrain.close();
        //==================================================================

//            //Read valid dataset to append to training dataset
//            //==================================================================
//            FileReader readerValid = new FileReader(pathDatasetValid);
//            BufferedReader rValid = new BufferedReader(readerValid);
//            line = null;
//            dataFound = 0;
//
//            while ((line = rValid.readLine()) != null) {
//
//                //Just read the file util do not find @DATA token
//                if (dataFound == 0) {
//
//                    Matcher mData = patternData.matcher(line);
//
//                    //See if reached @DATA token
//                    if (mData.find()) {
//                        dataFound = 1;
//                    }
//                } //Token @DATA was found
//                else {
//                    String[] vetLine = line.split(",");
//                    datasetValid.add(vetLine);
//                }
//            }
//            rValid.close();
//            readerValid.close();
//
//            datasetTrain.addAll(datasetValid);
//            datasetValid.clear();
//            //==================================================================
        //System.out.println();
    } catch (IOException ioe) {
        ioe.printStackTrace();
    }

    //Set means and modes to replace in case of missing attribute values
    setMeansAndModesTrain();
}

/* ===========================================================
     * Load the test data
     * =========================================================== */
public void readTestData(String pathDatasetTest) {
    datasetTest = new ArrayList<String[]>();
//        String pathDatasetTest = Paths.getPathDatasets() + Parameters.getFileDatasetTest();
    Pattern patternData = Pattern.compile(tokenData, Pattern.CASE_INSENSITIVE);

    try {
        FileReader readerTest = new FileReader(pathDatasetTest);
        BufferedReader rTest = new BufferedReader(readerTest);
        String line = null;
        int dataFound = 0;

        while ((line = rTest.readLine()) != null) {

            //Just read the file util do not find @DATA token
            if (dataFound == 0) {
                Matcher mData = patternData.matcher(line);

                //See if reached @DATA token
                if (mData.find()) {
                    dataFound = 1;
                }
            } //Token @DATA was found
            else {
                String[] vetLine = line.split(",");
                datasetTest.add(vetLine);
            }
        }
        rTest.close();
        readerTest.close();
    } catch (IOException ioe) {
        ioe.printStackTrace();
    }

    //Set means and modes to replace in case of missing attribute values
    //setMeansAndModesTest();
}

/* ===========================================================
     * Get the means and modes to substitute missing values
     * =========================================================== */
private void setMeansAndModesTrain() {
    int numAttributes = infoAttributes.size();

    //Train dataset
    //---------------------------------------------------
    int numExamples = datasetTrain.size();
    double[] replaceValues = new double[numAttributes];
    ArrayList<String> catReplaceValues = new ArrayList<String>();
    int numCategoricAttribute = 0;

    for (int i = 0; i < numAttributes; i++) {
        if (infoAttributes.get(i) == 1) { //numeric attribute
            //System.out.println("Num Attribute = " + i);
            double sum = 0;
            double mean = 0;
            int numNotMissing = 0;

            for (int j = 0; j < numExamples; j++) {
                //System.out.println("\tNum Example = " + j);

                String value = datasetTrain.get(j)[i];
                if (value.equals(tokenMissingValue) == false) {
                    sum += Double.parseDouble(value);
                    numNotMissing++;
                }
            }

            mean = sum / numNotMissing;
            replaceValues[i] = mean;

        } else { //categoric attibute
            replaceValues[i] = numCategoricAttribute;

            numCategoricAttribute++;

            ArrayList<String> domainValues = new ArrayList<String>();
            ArrayList<String> allValues = new ArrayList<String>();

            for (int j = 0; j < numExamples; j++) {
                String value = datasetTrain.get(j)[i];

                if (value.equals(tokenMissingValue) == false) {
                    domainValues.add(value);
                    allValues.add(value);
                }
            }

            //Eliminate duplicated classes
            HashSet hs = new HashSet();
            hs.addAll(domainValues);
            domainValues.clear();
            domainValues.addAll(hs);

            //Get the value that appears most often in the set
            String mostFrequent = getMostFrequentValue(domainValues, allValues);

            catReplaceValues.add(mostFrequent);
        }
    }

    meanValues.add(replaceValues);
    modeValues.add(catReplaceValues);
    //---------------------------------------------------

    //Valid dataset
    //---------------------------------------------------
    /*numExamples = datasetValid.size();
         replaceValues = new double[numAttributes];
         catReplaceValues = new ArrayList<String>();
         numCategoricAttribute = 0;
        
         for (int i = 0; i < numAttributes; i++) {
         if (infoAttributes.get(i) == 1) { //numeric attribute
        
         double sum = 0;
         double mean = 0;
         int numNotMissing = 0;
        
         for (int j = 0; j < numExamples; j++) {
        
         String value = datasetValid.get(j)[i];
         if (value.equals(tokenMissingValue) == false) {
         sum += Double.parseDouble(value);
         numNotMissing++;
         }
         }
        
         mean = sum / numNotMissing;
         replaceValues[i] = mean;
         } else { //categoric attibute
        
         replaceValues[i] = numCategoricAttribute;
        
         numCategoricAttribute++;
        
         ArrayList<String> domainValues = new ArrayList<String>();
         ArrayList<String> allValues = new ArrayList<String>();
        
         for (int j = 0; j < numExamples; j++) {
        
         String value = datasetValid.get(j)[i];
         if (value.equals(tokenMissingValue) == false) {
         domainValues.add(value);
         allValues.add(value);
         }
         }
        
         //Eliminate duplicated classes
         HashSet hs = new HashSet();
         hs.addAll(domainValues);
         domainValues.clear();
         domainValues.addAll(hs);
        
         //Get the value that appears most often in the set
         String mostFrequent = getMostFrequentValue(domainValues, allValues);
        
         catReplaceValues.add(mostFrequent);
         }
         }
        
         meanValues.add(replaceValues);
         modeValues.add(catReplaceValues);*/
}

private void setMeansAndModesTest() {
    int numAttributes = infoAttributes.size();

    int numExamples = datasetTest.size();
    double[] replaceValues = new double[numAttributes];
    ArrayList<String> catReplaceValues = new ArrayList<String>();
    int numCategoricAttribute = 0;

    for (int i = 0; i < numAttributes; i++) {
        if (infoAttributes.get(i) == 1) { //numeric attribute
            double sum = 0;
            double mean = 0;
            int numNotMissing = 0;

            for (int j = 0; j < numExamples; j++) {
                String value = datasetTest.get(j)[i];

                if (value.equals(tokenMissingValue) == false) {
                    sum += Double.parseDouble(value);
                    numNotMissing++;
                }
            }

            mean = sum / numNotMissing;
            replaceValues[i] = mean;

        } else { //categoric attibute
            replaceValues[i] = numCategoricAttribute;

            numCategoricAttribute++;

            ArrayList<String> domainValues = new ArrayList<String>();
            ArrayList<String> allValues = new ArrayList<String>();

            for (int j = 0; j < numExamples; j++) {
                String value = datasetTest.get(j)[i];

                if (value.equals(tokenMissingValue) == false) {
                    domainValues.add(value);
                    allValues.add(value);
                }
            }

            //Eliminate duplicated classes
            HashSet hs = new HashSet();
            hs.addAll(domainValues);
            domainValues.clear();
            domainValues.addAll(hs);

            //Get the value that appears most often in the set
            String mostFrequent = getMostFrequentValue(domainValues, allValues);

            catReplaceValues.add(mostFrequent);
        }
    }

    meanValues.add(replaceValues);
    modeValues.add(catReplaceValues);
}

/* ===========================================================
     * Get the most frequent categoric value of a given attribute
     * =========================================================== */
private String getMostFrequentValue(ArrayList<String> domainValues, ArrayList<String> allValues) {
    String mostFrequent = "";

    int[] counts = new int[domainValues.size()];

    //Count the elements
    for (int i = 0; i < domainValues.size(); i++) {
        for (int j = 0; j < allValues.size(); j++) {
            String domainValue = domainValues.get(i);
            String value = allValues.get(j);

            if (domainValue.equals(value)) {
                counts[i]++;
            }
        }
    }

    //Get the index of maximum
    int maxIndex = getMaxIndex(counts);
    mostFrequent = domainValues.get(maxIndex);

    return mostFrequent;
}

/* ===========================================================
     * Get the index of the maximum element in the array
     * =========================================================== */
private int getMaxIndex(int[] counts) {
    int maxIndex = 0;
    int max = counts[0];

    for (int i = 1; i < counts.length; i++) {
        if (counts[i] >= max) {
            max = counts[i];
            maxIndex = i;
        }
    }

    return maxIndex;
}

/* ===========================================================
     * Get the possible values for a categoric attribute
     * =========================================================== */
private String getCategoricValues(String line) {
    String[] values = line.split("\\{");
    String[] values2 = values[1].split("}");

    return values2[0];
}

/* ===========================================================
     * Get the name of the attribute
     * =========================================================== */
public String getAttributeName(String line, int typeAttribute) {
    String attributeName = "";

    if (typeAttribute == 1) {
        //Numeric attribute
        String[] vectorLine = line.split("(?i)ATTRIBUTE");
        String[] vectorLine2 = vectorLine[1].split("(?i)numeric");
        attributeName = vectorLine2[0].trim();

    } else {
        //Categoric attribute
        String[] vectorLine = line.split("(?i)ATTRIBUTE");
        String[] vectorLine2 = vectorLine[1].split("\\{");
        attributeName = vectorLine2[0].trim();
    }

    return attributeName;
}

/* ===========================================================
     * Map the possible values that can be put in a categoric clausule
     * These include combinations of values for the "in" operator
     * =========================================================== */
private ArrayList<String[]> mapPossibleCategoricValues(String stringValues) {
    possibleCombinations = new ArrayList<String[]>();

    String[] possibleValues = stringValues.split(",");
    String[] aux = new String[possibleValues.length + 1];

    aux[0] = "\\0";
    getCombinations(possibleValues, aux, 0);

    return possibleCombinations;
}

/* ===========================================================
     * Get all possible combinations of the categorical attributes
     * happening together
     * =========================================================== */
private void getCombinations(String[] possibleValues, String[] actual, int pos) {
    int size = getSize(actual);

    for (int i = pos; i < possibleValues.length; i++) {
        actual[size] = possibleValues[i];
        actual[size + 1] = "\\0";

        String[] newActual = new String[size + 1];
        System.arraycopy(actual, 0, newActual, 0, size + 1);

        possibleCombinations.add(newActual);
        getCombinations(possibleValues, actual, i + 1);
    }
}


/* ===========================================================
     * Auxiliary to get the size of a string.
     * The character '\\0' represents the end of the string
     * =========================================================== */
private int getSize(String[] actual) {
    int size = 0;

    for (int i = 0; i < actual.length; i++) {
        if ("\\0".equals(actual[i])) {
            size = i;
            break;
        }
    }

    return size;
}

/* ===========================================================
     * Verify if a given attribute value is present in a given example
     * =========================================================== */
public static boolean verifyAttributeValue(String attributeValue, String[] example) {
    boolean present = false;

    for (int i = 0; i < example.length; i++) {
        if (attributeValue.equals(example[i]) == true) {
            present = true;
            break;
        }
    }

    return present;
}

/* ===========================================================
     * Remove examples from the dataset, given their indexes
     * =========================================================== */
public static void removeTrainExamples(ArrayList<Integer> indexesCoveredExamples) {
    Collections.sort(indexesCoveredExamples, Collections.reverseOrder());

    for (int i = 0; i < indexesCoveredExamples.size(); i++) {
        datasetTrain.remove((int) indexesCoveredExamples.get(i));
    }
}

/* ===========================================================
     * Free memory of dataset test
     * =========================================================== */
public void freeTestDataset() {
    datasetTest.clear();
}

public static ArrayList<ArrayList<String[]>> getCategoricMapping() {
    return categoricMapping;
}

public static ArrayList<Integer> getCategoricAttributes() {
    return categoricAttributes;
}

public static ArrayList<String> getCategoricValues() {
    return categoricValues;
}

public static ArrayList<String[]> getDatasetTest() {
    return datasetTest;
}

public static ArrayList<String[]> getDatasetTrain() {
    return datasetTrain;
}

public static ArrayList<String[]> getDatasetValid() {
    return datasetValid;
}

public static ArrayList<Integer> getInfoAttributes() {
    return infoAttributes;
}

public static ArrayList<double[]> getMeanValues() {
    return meanValues;
}

public static ArrayList<ArrayList<String>> getModeValues() {
    return modeValues;
}

public static String getTokenMissingValue() {
    return tokenMissingValue;
}

public static String getTokenRootClass() {
    return tokenRootClass;
}

public static ArrayList<String> getNamesAttributes() {
    return namesAttributes;
}
}
