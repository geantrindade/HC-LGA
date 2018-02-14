package br.ufscar.hclga.classes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Population {

//Will store the generated population
private ArrayList<double[]> population;
//Will store the active terms for each rule
private ArrayList<ArrayList<Integer>> activeTerms;
//Use the train dataset to produce the initial individual
private ArrayList<String[]> datasetTrain;
//private HashMap<String, ArrayList<String[]>> classesAndInstances = new HashMap<>();
//private HashMap<String, ArrayList<String[]>> prevClassesAndInstances = new HashMap<>();

public Population(boolean multiObj) {
    if (multiObj == false) {
        ArrayList<Integer> infoAttributes = Datasets.getInfoAttributes();
        int numberAttributes = infoAttributes.size();
        population = new ArrayList<double[]>();
        activeTerms = new ArrayList<ArrayList<Integer>>();
        datasetTrain = Datasets.getDatasetTrain();

        /* ========================================================================================
         * This code below initiates each individual using seeding, where an example
         * is randomly chosen from the dataset and transformed in a rule.
         * With this, it is guaranteed that each rule covers at least one example
         * ========================================================================================
         */
        int posExample = 0;
        int pos = 0;

        Random generator = new Random();

//    for (int i = 0; i < datasetTrain.size(); i++) {
//        String[] instance = datasetTrain.get(i);
//        String key = instance[instance.length - 1];
//
//        if (!classesAndInstances.containsKey(key)) {
//            ArrayList<String[]> arrayL = new ArrayList<>();
//            arrayL.add(instance);
//            classesAndInstances.put(key, arrayL);
//            prevClassesAndInstances.put(key, new ArrayList<String[]>());
//
//        } else {
//            ArrayList<String[]> arrayL = classesAndInstances.get(key);
//            arrayL.add(instance);
//            classesAndInstances.put(key, arrayL);
//        }
//    }
//    int posClass = 0;
        while (population.size() < Parameters.getNumberInitialRules()) {
            ArrayList<Integer> activeTermsRule = new ArrayList<Integer>();
            double[] individual = new double[numberAttributes * 4];
            pos = 0;

//        boolean check = true;
//        //Randomly select an example from the dataset.
            posExample = generator.nextInt(datasetTrain.size());
            String[] example = datasetTrain.get(posExample);

//        String[] example = {};
//        int attempts = 0;
//
//        //iterates until a rule be generate for a specific class with no seed instance repeat
//        while (check) {
//            if (posClass == classesAndInstances.size()) {
//                posClass = 0;
//            }
//
//            String key = (String) classesAndInstances.keySet().toArray()[posClass];
//            ArrayList<String[]> possibleInstances = classesAndInstances.get(key);
//
//            //verifies if all instances from a class were already used
//            if (prevClassesAndInstances.get(key).size() < possibleInstances.size()) {
//                attempts = 0;
//                posExample = generator.nextInt(possibleInstances.size());
//
//                while (prevClassesAndInstances.get(key).contains(posExample)) {
//                    posExample = generator.nextInt(possibleInstances.size());
//                }
//
//                example = possibleInstances.get(posExample);
//                ArrayList<String[]> getArray = prevClassesAndInstances.get(key);
//                getArray.add(example);
//                prevClassesAndInstances.put(key, getArray);
//                check = false;
//            } else {
//                attempts++;
//            }
//
//            //reset records about instances once used, because they may not have been covered previously
//            if (attempts == 100) {
//                ArrayList<String[]> arrayAux = prevClassesAndInstances.get(key);
//                arrayAux.clear();
//                prevClassesAndInstances.put(key, arrayAux);
//
//            }
//            posClass++;
//        }
            //No empty rules will be created
            int emptyRule = 0;

            for (int j = 0; j < numberAttributes; j++) {

                //Randomly set if the clausule will be used or not
                individual[pos] = Operators.getInitialFlagValue(Parameters.getProbabilityUseClausule());
                emptyRule += individual[pos];
                //System.out.print(individual[pos] + " ");

                //Store the active terms
                if (individual[pos] == 1) {
                    activeTermsRule.add(pos);
                }

                if (infoAttributes.get(j) == 1) { //Numeric attribute

                    //Randomly select an numeric operator
                    individual[pos + 1] = Operators.getNumericOperator();
                    //Only generate relational rules
                    //individual[pos + 1] = 3;

                    //If the attribute value is missing, put the mean value
                    double attributeValue = 0.0;
                    if (Datasets.getTokenMissingValue().equals(example[j]) == true) {
                        attributeValue = Datasets.getMeanValues().get(0)[j];
                    } else {
                        attributeValue = Double.parseDouble(example[j]);
                    }

                    //Choose values for the clausule depending on the operator chosen
                    if (individual[pos + 1] == 0) { //attribute <= value
                        individual[pos + 2] = 0.0;
                        individual[pos + 3] = attributeValue;

                    } else if (individual[pos + 1] == 1) { //attribute >= value
                        individual[pos + 2] = attributeValue;
                        individual[pos + 3] = 0.0;

                    } else if (individual[pos + 1] == 2) { //value <= attribute <= value
                        double exampleValue = 0.0;
                        exampleValue = attributeValue;

                        double[] lowerHigherValue = getLowerHigherValue(exampleValue);
                        individual[pos + 2] = lowerHigherValue[0];
                        individual[pos + 3] = lowerHigherValue[1];

                    } else if (individual[pos + 1] == 3) { //attribute <= attribute
                        individual[pos + 2] = getPosAttributeLowerValue(posExample, j, attributeValue);
                        individual[pos + 3] = j;

                        //If the attribute value is an extreme, we use a normal rule
                        if (individual[pos + 2] == -1) {

                            if (generator.nextDouble() < 0.5) {//attribute <= value
                                individual[pos + 1] = 0.0;
                                individual[pos + 2] = 0.0;
                                individual[pos + 3] = attributeValue;
                            } else {
                                individual[pos + 1] = 1.0;
                                individual[pos + 2] = attributeValue;
                                individual[pos + 3] = 0.0;
                            }
                        }
                    }

                } else { //Categoric attribute

                    //Randomly select a categoric operator
                    individual[pos + 1] = Operators.getCategoricOperator();
                    int posAttribute = Datasets.getCategoricAttributes().indexOf(j);

                    String attributeValue = "";
                    //If the attribute value is missing, put the mode value
                    if (Datasets.getTokenMissingValue().equals(example[j]) == true) {
                        int posCatValue = (int) Datasets.getMeanValues().get(0)[j];
                        attributeValue = Datasets.getModeValues().get(0).get(posCatValue);
                    } else {
                        attributeValue = example[j];
                    }

                    if (individual[pos + 1] == 0) { // attribute = value

                        individual[pos + 2] = getCategoricToNumericValue(posAttribute, attributeValue);
                        individual[pos + 3] = 0.0;

                    } else if (individual[pos + 1] == 1) {// attribute != value

                        for (int k = 0; k < Datasets.getCategoricMapping().get(posAttribute).size(); k++) {
                            if (Datasets.getCategoricMapping().get(posAttribute).get(k).length == 1) {
                                if (Datasets.getCategoricMapping().get(posAttribute).get(k)[0].equals(attributeValue) == false) {
                                    individual[pos + 2] = (double) k;
                                    individual[pos + 3] = 0.0;
                                    break;
                                }
                            }
                        }
                    } else if (individual[pos + 1] == 2) {// attribute in values

                        int found = 0;
                        //ArrayList<String[]> teste = Datasets.getCategoricMapping().get(posAttribute);
                        for (int k = 0; k < Datasets.getCategoricMapping().get(posAttribute).size(); k++) {
                            if (found == 1) {
                                break;
                            }
                            if (Datasets.getCategoricMapping().get(posAttribute).get(k).length > 1) {
                                for (int l = 0; l < Datasets.getCategoricMapping().get(posAttribute).get(k).length; l++) {
                                    //String[] teste = Datasets.getCategoricMapping().get(posAttribute).get(k);
                                    if (Datasets.getCategoricMapping().get(posAttribute).get(k)[l].equals(attributeValue) == true) {
                                        individual[pos + 2] = (double) k;
                                        individual[pos + 3] = 0.0;
                                        found = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                pos += 4;
            }

            //System.out.println();
            //Add generated individual to the population
            if (emptyRule > 0) {
                population.add(individual);
                activeTerms.add(activeTermsRule);
            }
        }
    }
}

/* ===========================================================
     * Given the indexes of an example and an attribute, searches
     * for another attribute that has a lower or equal value and
     * returns its index
     * =========================================================== */
public int getPosAttributeLowerValue(int posExample, int posAttribute, double attributeValue) {

    int pos = -1;
    int numberAttributes = Datasets.getInfoAttributes().size();

    for (int i = 0; i < numberAttributes; i++) {
        if (Datasets.getInfoAttributes().get(i) == 1 && i != posAttribute) {
            double searchValue = 0;

            if (Datasets.getTokenMissingValue().equals(getDatasetTrain().get(posExample)[i]) == true) {
                searchValue = Datasets.getMeanValues().get(0)[i];
            } else {
                searchValue = Double.parseDouble(getDatasetTrain().get(posExample)[i]);
            }

            if (searchValue <= attributeValue) {
                pos = i;
                break;
            }
        }
    }
    /*if (pos > 80) {
            System.out.println();
        }*/

    return pos;
}

/* ===========================================================
     * Given an categoric attribute value and the number of the
     * attribute, returns the corresponding numeric value of this
     * categoric value
     * =========================================================== */
public static double getCategoricToNumericValue(int posAttribute, String attributeValue) {

    double numericValue = 0.0;

    for (int k = 0; k < Datasets.getCategoricMapping().get(posAttribute).size(); k++) {
        if (Datasets.getCategoricMapping().get(posAttribute).get(k).length == 1) {
            if (Datasets.getCategoricMapping().get(posAttribute).get(k)[0].equals(attributeValue) == true) {
                numericValue = (double) k;
                break;
            }
        }
    }

    return numericValue;
}

/* ===========================================================
     * Get a value for clausules (x in {value1, value2, ..., valuen})
     * =========================================================== */
private double getSetOfCategoricValues(ArrayList<Integer> categoricValues) {

    double value = 0.0;

    Random generator = new Random();
    int pos = generator.nextInt(categoricValues.size());
    value = categoricValues.get(pos);

    return value;

}

/* ===========================================================
     * Get a value for clausules (x = value) or (x != value)
     * =========================================================== */
private double getCategoricValue(ArrayList<Integer> categoricValues) {

    double value = 0.0;

    Random generator = new Random();
    int pos = generator.nextInt(categoricValues.size());
    value = categoricValues.get(pos);

    return value;

}

/* ===========================================================
     * Get a value for clausules (x < value) or (x <= value)
     * or (x < value) or (x <= value)
     * =========================================================== */
private double getNumericValue(int numAttribute) {

    double value = 0.0;

    Random generator = new Random();
    int numInstance = generator.nextInt(getDatasetTrain().size());
    String stringValue = getDatasetTrain().get(numInstance)[numAttribute];

    while (stringValue.equals(Datasets.getTokenMissingValue()) == true) {
        numInstance = generator.nextInt(getDatasetTrain().size());
        stringValue = getDatasetTrain().get(numInstance)[numAttribute];
    }

    value = Double.parseDouble(stringValue);

    return value;
}

/* ===========================================================
     * Get a value for clausules (value1 < x <= value2)
     * or (value1 <= x < value2)
     * =========================================================== */
 /*private double[] getLowerHigherValue(int numAttribute) {
    
     double[] lowerHigherValue = {0.0, 0.0};
    
     while (lowerHigherValue[0] >= lowerHigherValue[1]) {
    
     Random generator = new Random();
    
     int numInstance0 = generator.nextInt(datasetTrain.size());
     String stringValue0 = datasetTrain.get(numInstance0)[numAttribute];
    
     int numInstance1 = generator.nextInt(datasetTrain.size());
     String stringValue1 = datasetTrain.get(numInstance1)[numAttribute];
    
     while (stringValue0.equals(Datasets.getTokenMissingValue()) == true || stringValue1.equals(Datasets.getTokenMissingValue()) == true) {
    
     numInstance0 = generator.nextInt(datasetTrain.size());
     stringValue0 = datasetTrain.get(numInstance0)[numAttribute];
    
     numInstance1 = generator.nextInt(datasetTrain.size());
     stringValue1 = datasetTrain.get(numInstance1)[numAttribute];
    
     }
    
     lowerHigherValue[0] = Double.parseDouble(stringValue0);
     lowerHigherValue[1] = Double.parseDouble(stringValue1);
    
     }
    
     return lowerHigherValue;
     }*/
private double[] getLowerHigherValue(double exampleValue) {

    double[] lowerHigherValue = {0.0, 0.0};

    while (lowerHigherValue[0] > exampleValue || lowerHigherValue[1] < exampleValue) {

        Random generator = new Random();

        if (exampleValue != 0) {

            lowerHigherValue[0] = (4 * exampleValue) * generator.nextDouble() - exampleValue;
            lowerHigherValue[1] = (4 * exampleValue) * generator.nextDouble() - exampleValue;
        } else {

            lowerHigherValue[0] = 2 * generator.nextDouble() - 1;
            lowerHigherValue[1] = 2 * generator.nextDouble() - 1;

        }

    }

    return lowerHigherValue;
}

/*
     private double[] getLowerHigherValue(int numAttribute, double exampleValue) {
    
     double[] lowerHigherValue = {0.0, 0.0};
    
     while (lowerHigherValue[0] > exampleValue || lowerHigherValue[1] < exampleValue) {
    
     Random generator = new Random();
    
     int numInstance0 = generator.nextInt(datasetTrain.size());
     String stringValue0 = datasetTrain.get(numInstance0)[numAttribute];
    
     int numInstance1 = generator.nextInt(datasetTrain.size());
     String stringValue1 = datasetTrain.get(numInstance1)[numAttribute];
    
     while (stringValue0.equals(Datasets.getTokenMissingValue()) == true || stringValue1.equals(Datasets.getTokenMissingValue()) == true) {
    
     numInstance0 = generator.nextInt(datasetTrain.size());
     stringValue0 = datasetTrain.get(numInstance0)[numAttribute];
    
     numInstance1 = generator.nextInt(datasetTrain.size());
     stringValue1 = datasetTrain.get(numInstance1)[numAttribute];
    
     }
    
     lowerHigherValue[0] = Double.parseDouble(stringValue0);
     lowerHigherValue[1] = Double.parseDouble(stringValue1);
    
     }
    
     return lowerHigherValue;
     }
     * 
 */
public ArrayList<double[]> getPopulation() {
    return population;
}

public ArrayList<ArrayList<Integer>> getActiveTerms() {
    return activeTerms;
}

    /**
     * @param population the population to set
     */
    public void setPopulation(ArrayList<double[]> population) {
        this.population = population;
    }

    /**
     * @param activeTerms the activeTerms to set
     */
    public void setActiveTerms(ArrayList<ArrayList<Integer>> activeTerms) {
        this.activeTerms = activeTerms;
    }

    /**
     * @return the datasetTrain
     */
    public ArrayList<String[]> getDatasetTrain() {
        return datasetTrain;
    }

    /**
     * @param datasetTrain the datasetTrain to set
     */
    public void setDatasetTrain(ArrayList<String[]> datasetTrain) {
        this.datasetTrain = datasetTrain;
    }

}
