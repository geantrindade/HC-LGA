package br.ufscar.hclga.classes;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Classes {

private static ArrayList<int[]> binaryClassesTrain;
private static ArrayList<int[]> binaryClassesValid;
private static ArrayList<int[]> binaryClassesTest;
private static ArrayList<ArrayList<String>> DAGrelationships;
private static String[] classes;
private static double[] weightingScheme;
private static double[] meanClassLabelVectorAllClasses;
private static ArrayList<ArrayList<Integer>> positionClassesLevel;
//Illegal GO classes
private static String[] illegalGOclasses = {"GO0003674", "GO0005575", "GO0008150"};
//private static String[] illegalGOclasses = {};

/* ===========================================================
* Set a vector with all classes of the DAG structure
* =========================================================== */
public static void setDAGclasses(String lineClasses, String tokenHierarchical) {
    DAGrelationships = new ArrayList<ArrayList<String>>();
    ArrayList<String> classesAux = new ArrayList<String>();
    String[] vetLine = lineClasses.split(tokenHierarchical);
    String[] classesAux2 = vetLine[1].split(",");

    //These arrays will store the class relationships
    //parents.get(pos) is the superclass of children.get(pos)
    ArrayList<String> parents = new ArrayList<String>();
    ArrayList<String> children = new ArrayList<String>();

    classesAux.add(classesAux2[0].split("/")[1].trim());
    parents.add(classesAux2[0].split("/")[0].trim());
    children.add(classesAux2[0].split("/")[1].trim());

    for (int i = 1; i < classesAux2.length; i++) {
        String[] classesAux3 = classesAux2[i].split("/");
        
        parents.add(classesAux3[0].trim());
        children.add(classesAux3[1].trim());

        if (!classesAux.contains(classesAux3[0]) && !Datasets.getTokenRootClass().equals(classesAux3[0])) {
            classesAux.add(classesAux3[0].trim());
        }

        if (!classesAux.contains(classesAux3[1])) {
            classesAux.add(classesAux3[1].trim());
        }
    }

    //Vector with all classes
    classes = new String[classesAux.size()];

    for (int i = 0; i < classesAux.size(); i++) {
        classes[i] = classesAux.get(i);
    }

    //Class relationships
    DAGrelationships.add(parents);
    DAGrelationships.add(children);

    //Weighting scheme for the classes
    weightingScheme = new double[classes.length];
    setWeightingScheme(weightingScheme);
}

/* ===========================================================
* Set a vector with all classes of the tree structure
* =========================================================== */
public static void setTreeClasses(String lineClasses, String tokenHierarchical) {
    String[] vetLine = lineClasses.split(tokenHierarchical);
    classes = vetLine[1].split(",");
    classes[0] = classes[0].trim();

    //Weighting scheme for the classes
    weightingScheme = new double[classes.length];
    setWeightingScheme(weightingScheme);

    //Set the position of the classes by level
    positionClassesLevel = setPositionClassesLevel();
}

/*===========================================================================
* Get position of classes by level
*===========================================================================*/
public static ArrayList<ArrayList<Integer>> setPositionClassesLevel() {
    ArrayList<ArrayList<Integer>> positionClassesLevels = new ArrayList<ArrayList<Integer>>();
    String rootClass = "";

    for (int i = 0; i < Parameters.getNumLevels(); i++) {
        ArrayList<Integer> positions = new ArrayList<Integer>();
        rootClass = rootClass.concat("[0-9]+");
        Pattern pattern = Pattern.compile("^" + rootClass + "$");

        for (int j = 0; j < classes.length; j++) {
            Matcher m = pattern.matcher(classes[j]);

            if (m.find()) {
                positions.add(j);
            }
        }

        rootClass = rootClass.concat("/");
        positionClassesLevels.add(positions);
    }

    return positionClassesLevels;
}

/* ===========================================================
* Set the mean class label vector considering all classes of the dataset
* =========================================================== */
public static void setMeanClassLabelVectorAll() {
    for (int i = 0; i < binaryClassesTrain.size(); i++) {
        for (int j = 0; j < classes.length; j++) {
            meanClassLabelVectorAllClasses[j] += binaryClassesTrain.get(i)[j];
        }
    }

    for (int i = 0; i < meanClassLabelVectorAllClasses.length; i++) {
        meanClassLabelVectorAllClasses[i] = meanClassLabelVectorAllClasses[i] / binaryClassesTrain.size();
    }
}

/* ===========================================================
* Set weights for the classes according to Vens et al, 2008
* Decision Trees for Hierarchical Multi-Label Classification
* Machine Learning 73(2):185-214
* =========================================================== */
private static void setWeightingScheme(double[] weightingScheme) {
    //The top level classes of the hierarchy receive a weigthing value of 0.75
    ArrayList<Integer> topLevelClassesPositions = getTopLevelClassesPositions();

    if (Parameters.getHierarchyType().equals("Tree")) {//Tree hierarchy
        for (int i = 0; i < topLevelClassesPositions.size(); i++) {
            weightingScheme[topLevelClassesPositions.get(i)] = 0.75;
            String rootClass = "^" + classes[topLevelClassesPositions.get(i)];
            ArrayList<Integer> superClassesPositions = new ArrayList<Integer>();
            superClassesPositions.add(topLevelClassesPositions.get(i));
            setTreeClassesWeights(weightingScheme, rootClass, superClassesPositions);
        }
        
    } else {//DAG hierarchy
        for (int i = 0; i < topLevelClassesPositions.size(); i++) {
            weightingScheme[topLevelClassesPositions.get(i)] = 0.75;
            String rootClass = classes[topLevelClassesPositions.get(i)];
            ArrayList<Integer> superClassesPositions = new ArrayList<Integer>();
            superClassesPositions.add(topLevelClassesPositions.get(i));
            setDAGClassesWeights(weightingScheme, rootClass);
        }
    }
}

/* ===========================================================
* Recursive method to set the weights of the DAG classes
* given a root class
* =========================================================== */
private static void setDAGClassesWeights(double[] weightingScheme, String rootClass) {
    //Get the closest children of the root class
    ArrayList<String> children = new ArrayList<String>();

    for (int i = 0; i < DAGrelationships.get(0).size(); i++) {
        if (DAGrelationships.get(0).get(i).equals(rootClass) == true) {
            children.add(DAGrelationships.get(1).get(i));
        }
    }

    //Get all parent classes of each children class, as a class
    //can have more than one superclass
    for (int i = 0; i < children.size(); i++) {
        //Gets the position of this child class in classes
        int posChild = 0;
        ArrayList<String> parents = new ArrayList<String>();
        HashSet hs = new HashSet();
        ArrayList<Integer> superClassesPositions = new ArrayList<Integer>();

        for (int j = 0; j < classes.length; j++) {
            if (children.get(i).equals(classes[j]) == true) {
                posChild = j;
            }
        }

        parents.add(children.get(i));
        getAllPossibleDAGclasses(parents, children.get(i));

        //Eliminate duplicated classes
        hs.addAll(parents);
        parents.clear();
        parents.addAll(hs);
        parents.remove(Datasets.getTokenRootClass());
        parents.remove(children.get(i));

        //Verify if the weight of any parent is not set.
        //If a parent weight is not set, we go to the next child class
        int notSet = 0;

        for (int j = 0; j < parents.size(); j++) {
            for (int k = 0; k < classes.length; k++) {
                if (parents.get(j).equals(classes[k]) == true) {
                    //Get positions of the superclasses
                    superClassesPositions.add(k);
                    if (weightingScheme[k] == 0.0) {
                        notSet = 1;
                        break;
                    }
                }
            }

            if (notSet == 1) {
                break;
            }
        }

        //If all parent weights are set, calculate the weight
        //and call the method using the actual child class
        if (notSet == 0) {
            int numParents = 0;
            double sum = 0.0;

            for (int j = 0; j < superClassesPositions.size(); j++) {
                numParents++;
                sum += weightingScheme[superClassesPositions.get(j)];
            }

            weightingScheme[posChild] = 0.75 * (sum / numParents);
            setDAGClassesWeights(weightingScheme, children.get(i));
        }
    }
}

/* ===========================================================
* Recursive method to set the weights of the Tree classes
* given a root class
* =========================================================== */
private static void setTreeClassesWeights(double[] weightingScheme, String rootClass, ArrayList<Integer> superClassesPositionsAux) {
    //Get the children of the root class in the next level
    rootClass = rootClass.concat("/[0-9]+");
    Pattern pattern = Pattern.compile(rootClass + "$");

    for (int i = 0; i < classes.length; i++) {
        Matcher m = pattern.matcher(classes[i]);

        if (m.find()) {
            int numParents = 0;
            double sum = 0.0;
            ArrayList<Integer> superClassesPositions = new ArrayList<Integer>();

            for (int j = 0; j < superClassesPositionsAux.size(); j++) {
                numParents++;
                sum += weightingScheme[superClassesPositionsAux.get(j)];
                superClassesPositions.add(superClassesPositionsAux.get(j));
            }

            weightingScheme[i] = 0.75 * (sum / numParents);
            superClassesPositions.add(i);
            setTreeClassesWeights(weightingScheme, "^" + classes[i], superClassesPositions);
        }
    }
}

/* ===========================================================
* Get the top level classes of the hierarchy
* =========================================================== */
public static ArrayList<Integer> getTopLevelClassesPositions() {
    ArrayList<Integer> topLevelClassesPositions = new ArrayList<Integer>();

    if (Parameters.getHierarchyType().equals("Tree")) {//Tree hierarchy
        for (int i = 0; i < classes.length; i++) {
            if (classes[i].contains("/") == false) {
                topLevelClassesPositions.add(i);
            }
        }
        
    } else {//DAG hierarchy
        for (int i = 0; i < DAGrelationships.get(0).size(); i++) {
            if (DAGrelationships.get(0).get(i).equals(Datasets.getTokenRootClass()) == true) {
                for (int j = 0; j < classes.length; j++) {
                    if (DAGrelationships.get(1).get(i).equals(classes[j]) == true) {
                        topLevelClassesPositions.add(j);
                        break;
                    }
                }
            }
        }
    }

    return topLevelClassesPositions;
}

/* ===========================================================
* Build the binary structure to store the dataset's classes
* =========================================================== */
public static void buildClassesStructureTrain() {
    int numberClasses = 0;
    numberClasses = classes.length;
    binaryClassesTrain = new ArrayList<int[]>();

    //Training dataset
    ArrayList<String[]> datasetTrain = Datasets.getDatasetTrain();

    for (int i = 0; i < datasetTrain.size(); i++) {
        int[] binaryVector = new int[numberClasses];
        String actualClasses = datasetTrain.get(i)[datasetTrain.get(i).length - 1];
        ArrayList<Integer> posClasses = getPosClasses(actualClasses);

        for (int j = 0; j < posClasses.size(); j++) {
            binaryVector[posClasses.get(j)] = 1;
        }

        binaryClassesTrain.add(binaryVector);
    }

    //Set the mean class label vector considering all classes
    meanClassLabelVectorAllClasses = new double[classes.length];
    setMeanClassLabelVectorAll();
}

public static void buildClassesStructureValid() {
    int numberClasses = 0;
    numberClasses = classes.length;
    binaryClassesValid = new ArrayList<int[]>();

    //Valid dataset
    ArrayList<String[]> datasetValid = Datasets.getDatasetValid();
    int[] binaryVector = new int[numberClasses];
    String actualClasses = "";
    ArrayList<Integer> posClasses;

    for (int i = 0; i < datasetValid.size(); i++) {
        actualClasses = datasetValid.get(i)[datasetValid.get(i).length - 1];
        posClasses = getPosClasses(actualClasses);

        for (int j = 0; j < posClasses.size(); j++) {
            binaryVector[posClasses.get(j)] = 1;
        }

        binaryClassesValid.add(binaryVector);
    }
}

public static void buildClassesStructureTest() {
    int numberClasses = 0;
    numberClasses = classes.length;
    binaryClassesTest = new ArrayList<int[]>();

    //Test dataset
    ArrayList<String[]> datasetTest = Datasets.getDatasetTest();

    for (int i = 0; i < datasetTest.size(); i++) {
        int[] binaryVector = new int[numberClasses];
        String actualClasses = datasetTest.get(i)[datasetTest.get(i).length - 1];
        ArrayList<Integer> posClasses = getPosClasses(actualClasses);

        for (int j = 0; j < posClasses.size(); j++) {
            binaryVector[posClasses.get(j)] = 1;
        }

        binaryClassesTest.add(binaryVector);
    }

    String directory = "Evaluation/";
    String file = "realTestClasses.txt";
    new File(directory).mkdirs();
    BufferedWriter writer;

    try {
        writer = new BufferedWriter(new FileWriter(directory + file));

        for (int i = 0; i < binaryClassesTest.size(); i++) {
            for (int j = 0; j < numberClasses; j++) {
                writer.write(binaryClassesTest.get(i)[j] + " ");
            }

            writer.newLine();
        }

        writer.flush();
        writer.close();

    } catch (IOException ex) {
        Logger.getLogger(Classes.class.getName()).log(Level.SEVERE, null, ex);
    }
}

/* ===========================================================
* Get the positions given classes in the binary vector
* =========================================================== */
public static ArrayList<Integer> getPosClasses(String actualClasses) {
    ArrayList<Integer> positions = new ArrayList<Integer>();
    String[] vectorClasses = actualClasses.split("@");
    ArrayList<String> allClasses = getAllPossibleClasses(vectorClasses);

    if ("Tree".equals(Parameters.getHierarchyType())) {
        for (int i = 0; i < allClasses.size(); i++) {
            for (int j = 0; j < classes.length; j++) {
                if (allClasses.get(i).equals(classes[j]) == true) {
                    positions.add(j);
                    break;
                }
            }
        }
        
    } else if ("DAG".equals(Parameters.getHierarchyType())) {
        allClasses.remove(Datasets.getTokenRootClass());

        for (int i = 0; i < allClasses.size(); i++) {
            for (int j = 0; j < classes.length; j++) {
                if (allClasses.get(i).equals(classes[j]) == true) {
                    positions.add(j);
                    break;
                }
            }
        }
    }

    return positions;
}

/* ===========================================================
* Given a vector with classes, return all classes in all levels
* =========================================================== */
public static ArrayList<String> getAllPossibleClasses(String[] vectorClasses) {
    ArrayList<String> allClasses = new ArrayList<String>();

    if ("Tree".equals(Parameters.getHierarchyType())) {
        String[] vetClasses;
        String aClass = "";

        for (int i = 0; i < vectorClasses.length; i++) {
            vetClasses = vectorClasses[i].split("/");

            for (int j = 0; j < vetClasses.length; j++) {
                aClass = aClass.concat(vetClasses[j]);
                allClasses.add(aClass);

                if (j < vetClasses.length - 1) {
                    aClass = aClass.concat("/");
                }
            }
        }
        
    } else if ("DAG".equals(Parameters.getHierarchyType())) {
        ArrayList<String> DAGclasses = new ArrayList<String>();

        for (int i = 0; i < vectorClasses.length; i++) {
            DAGclasses.add(vectorClasses[i]);
            getAllPossibleDAGclasses(DAGclasses, vectorClasses[i]);
            allClasses.addAll(DAGclasses);
        }
    }

    //Eliminate duplicated classes
    HashSet hs = new HashSet();
    hs.addAll(allClasses);
    allClasses.clear();
    allClasses.addAll(hs);

    return allClasses;
}

/* ===========================================================
* Recursive function to get parent classes of a given DAG class
* =========================================================== */
private static void getAllPossibleDAGclasses(ArrayList<String> DAGclasses, String DAGclass) {
    ArrayList<String> parents = DAGrelationships.get(0);
    ArrayList<String> children = DAGrelationships.get(1);

    for (int i = 0; i < children.size(); i++) {
        if (children.get(i).equals(DAGclass) == true) {
            DAGclasses.add(parents.get(i));
            getAllPossibleDAGclasses(DAGclasses, parents.get(i));
        }
    }
}

public static String[] getClasses() {
    return classes;
}

public static ArrayList<int[]> getBinaryClassesTest() {
    return binaryClassesTest;
}

public static ArrayList<int[]> getBinaryClassesTrain() {
    return binaryClassesTrain;
}

public static ArrayList<int[]> getBinaryClassesValid() {
    return binaryClassesValid;
}

public static double[] getWeightingScheme() {
    return weightingScheme;
}

public static double[] getMeanClassLabelVectorAllClasses() {
    return meanClassLabelVectorAllClasses;
}

public static String[] getIllegalGOclasses() {
    return illegalGOclasses;
}

public static ArrayList<ArrayList<Integer>> getPositionClassesLevel() {
    return positionClassesLevel;
}

}
