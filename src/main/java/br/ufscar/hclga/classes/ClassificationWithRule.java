package br.ufscar.hclga.classes;

import static br.ufscar.hclga.classes.Evaluation.applyThresholds;
import static br.ufscar.hclga.classes.Evaluation.getDataInterpolation;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class ClassificationWithRule {

private String nameDataset = leConfig()[0];
private int specificFold = Integer.parseInt(leConfig()[1]);
private int folds = Integer.parseInt(leConfig()[2]);
private String pathTrainDataset = leConfig()[3];
private String pathTestDataset = leConfig()[4];
private String pathRules = leConfig()[5];
private int numberOfRules;
private String pathToSavePredictions = leConfig()[6];

private ArrayList<String[]> datasetTrain;
private ArrayList<String[]> datasetTest;
private ArrayList<int[]> binaryClassesTrain;
private ArrayList<int[]> binaryClassesTest;
private String[] classes;
private double[] weightingScheme;
private double[] meanClassLabelVectorAllClassesTrain;
private double[] meanClassLabelVectorAllClassesTest;
private static ArrayList<ArrayList<Integer>> positionClassesLevel;
private int numberOfClasses;
private String[] attributes = {"AA", "AT", "AC", "AG", "TT", "TA", "TC", "TG", "CC", "CA", "CG", "CT", "GG", "GA", "GT", "GC", "AAA", "AAT", "AAC", "AAG", "ATA", "ACA", "AGA", "ATT", "ATC", "ATG",
    "ACC", "ACT", "ACG", "AGG", "AGC", "AGT", "CCC", "CCT", "CCA", "CCG", "CTC", "CAC", "CGC", "CTT", "CTA", "CTG", "CGG", "CGA", "CGT", "CAA", "CAT", "CAG", "GGG", "GGT", "GGA", "GGC", "GTG", "GAG", "GCG", "GTT", "GTA", "GTC",
    "GCC", "GCA", "GCT", "GAA", "GAT", "GAC", "TTT", "TTA", "TTC", "TTG", "TAT", "TCT", "TGT", "TGG", "TGA", "TGC", "TAA", "TAC", "TAG", "TCC", "TCA", "TCG", "AAAA", "AAAT", "AAAC", "AAAG", "AATA", "AACA", "AAGA", "AATT", "AATC", "AATG",
    "AACC", "AACT", "AACG", "AAGG", "AAGC", "AAGT", "ACCC", "ACCT", "ACCA", "ACCG", "ACTC", "ACAC", "ACGC", "ACTT", "ACTA", "ACTG", "ACGG", "ACGA", "ACGT", "ACAA", "ACAT", "ACAG", "AGGG", "AGGT", "AGGA", "AGGC", "AGTG", "AGAG", "AGCG", "AGTT", "AGTA", "AGTC",
    "AGCC", "AGCA", "AGCT", "AGAA", "AGAT", "AGAC", "ATTT", "ATTA", "ATTC", "ATTG", "ATAT", "ATCT", "ATGT", "ATGG", "ATGA", "ATGC", "ATAA", "ATAC", "ATAG", "ATCC", "ATCA", "ATCG", "CAAA", "CAAT", "CAAC", "CAAG", "CATA", "CACA", "CAGA", "CATT", "CATC", "CATG",
    "CACC", "CACT", "CACG", "CAGG", "CAGC", "CAGT", "CCCC", "CCCT", "CCCA", "CCCG", "CCTC", "CCAC", "CCGC", "CCTT", "CCTA", "CCTG", "CCGG", "CCGA", "CCGT", "CCAA", "CCAT", "CCAG", "CGGG", "CGGT", "CGGA", "CGGC", "CGTG", "CGAG", "CGCG", "CGTT", "CGTA", "CGTC",
    "CGCC", "CGCA", "CGCT", "CGAA", "CGAT", "CGAC", "CTTT", "CTTA", "CTTC", "CTTG", "CTAT", "CTCT", "CTGT", "CTGG", "CTGA", "CTGC", "CTAA", "CTAC", "CTAG", "CTCC", "CTCA", "CTCG", "GAAA", "GAAT", "GAAC", "GAAG", "GATA", "GACA", "GAGA", "GATT", "GATC", "GATG",
    "GACC", "GACT", "GACG", "GAGG", "GAGC", "GAGT", "GCCC", "GCCT", "GCCA", "GCCG", "GCTC", "GCAC", "GCGC", "GCTT", "GCTA", "GCTG", "GCGG", "GCGA", "GCGT", "GCAA", "GCAT", "GCAG", "GGGG", "GGGT", "GGGA", "GGGC", "GGTG", "GGAG", "GGCG", "GGTT", "GGTA", "GGTC",
    "GGCC", "GGCA", "GGCT", "GGAA", "GGAT", "GGAC", "GTTT", "GTTA", "GTTC", "GTTG", "GTAT", "GTCT", "GTGT", "GTGG", "GTGA", "GTGC", "GTAA", "GTAC", "GTAG", "GTCC", "GTCA", "GTCG", "TAAA", "TAAT", "TAAC", "TAAG", "TATA", "TACA", "TAGA", "TATT", "TATC", "TATG",
    "TACC", "TACT", "TACG", "TAGG", "TAGC", "TAGT", "TCCC", "TCCT", "TCCA", "TCCG", "TCTC", "TCAC", "TCGC", "TCTT", "TCTA", "TCTG", "TCGG", "TCGA", "TCGT", "TCAA", "TCAT", "TCAG", "TGGG", "TGGT", "TGGA", "TGGC", "TGTG", "TGAG", "TGCG", "TGTT", "TGTA", "TGTC",
    "TGCC", "TGCA", "TGCT", "TGAA", "TGAT", "TGAC", "TTTT", "TTTA", "TTTC", "TTTG", "TTAT", "TTCT", "TTGT", "TTGG", "TTGA", "TTGC", "TTAA", "TTAC", "TTAG", "TTCC", "TTCA", "TTCG"};

public ClassificationWithRule() {
    if (nameDataset.equalsIgnoreCase("Mips")) {
        numberOfClasses = 14;
    } else if (nameDataset.equalsIgnoreCase("Repbase")) {
        numberOfClasses = 31;
    } else {
        numberOfClasses = -1;
    }
//
//    setMeanClassLabelVectorAll();
}

public ArrayList<String> readRulesFile(String rulesFile) {
    ArrayList<String> result = new ArrayList<>();

    try {
        FileReader reader = new FileReader(rulesFile);
        BufferedReader buffReader = new BufferedReader(reader);

        String line = null;
        while ((line = buffReader.readLine()) != null) {
            if (!line.isEmpty()) {
                result.add(line);
            }
        }

        buffReader.close();
        reader.close();

    } catch (IOException ioe) {
        ioe.printStackTrace();
    }

    return result;
}

public ArrayList<String> readRulesFileClus(String rulesFile) {
    ArrayList<String> result = new ArrayList<>();
    String wholeRule = "";
    try {
        FileReader reader = new FileReader(rulesFile);
        BufferedReader buffReader = new BufferedReader(reader);

        String line = null;
        while ((line = buffReader.readLine()) != null) {
            if (!line.isEmpty()) {
                wholeRule += line;

                if (line.contains("THEN")) {
                    result.add(wholeRule);
                    wholeRule = "";
                }
            }
        }

        buffReader.close();
        reader.close();

    } catch (IOException ioe) {
        ioe.printStackTrace();
    }

    for (int i = 0; i < result.size(); i++) {
        String aux = result.get(i);
        aux = aux.replaceAll("=", "");
        aux = aux.replace(":", "=");
        aux = aux.replace("IF", "");
        aux = aux.replace("[", "");
        int index = aux.indexOf("]");
        aux = aux.substring(0, index);
        aux = aux.trim();
        result.remove(i);
        result.add(i, aux);
    }

    return result;
}

public void readTestData(String testDatasetFile) {
    datasetTest = new ArrayList<String[]>();
    int numAttribute = -1;

    Pattern patternData = Pattern.compile("@Data", Pattern.CASE_INSENSITIVE);
    Pattern patternNumeric = Pattern.compile("numeric", Pattern.CASE_INSENSITIVE);
    Pattern patternAttribute = Pattern.compile("@ATTRIBUTE", Pattern.CASE_INSENSITIVE);
    Pattern patternHierarchical = Pattern.compile("hierarchical", Pattern.CASE_INSENSITIVE);

    try {
        FileReader readerTest = new FileReader(testDatasetFile);
        BufferedReader rTest = new BufferedReader(readerTest);
        String line = null;
        int dataFound = 0;

        while ((line = rTest.readLine()) != null) {

            //Just read the file util do not find @DATA token
            if (dataFound == 0) {
                Matcher mAttribute = patternAttribute.matcher(line);
                Matcher mData = patternData.matcher(line);

                //See if reached @DATA token
                if (mData.find()) {
                    dataFound = 1;

                } else if (mAttribute.find()) {
                    numAttribute++;

                    //If so, check if attribute is hierarchical, numeric or categoric
                    Matcher mNumeric = patternNumeric.matcher(line);
                    Matcher mHierarchical = patternHierarchical.matcher(line);

                    if (mHierarchical.find()) {
                        //Build structure to store the classes
                        setTreeClasses(line, "hierarchical");
                    }
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
}

public void readTrainData(String trainDatasetFile) {
    datasetTrain = new ArrayList<String[]>();
    int numAttribute = -1;

    Pattern patternData = Pattern.compile("@Data", Pattern.CASE_INSENSITIVE);
    Pattern patternNumeric = Pattern.compile("numeric", Pattern.CASE_INSENSITIVE);
    Pattern patternAttribute = Pattern.compile("@ATTRIBUTE", Pattern.CASE_INSENSITIVE);
    Pattern patternHierarchical = Pattern.compile("hierarchical", Pattern.CASE_INSENSITIVE);

    try {
        FileReader readerTrain = new FileReader(trainDatasetFile);
        BufferedReader rTrain = new BufferedReader(readerTrain);
        String line = null;
        int dataFound = 0;

        while ((line = rTrain.readLine()) != null) {

            //Just read the file util do not find @DATA token
            if (dataFound == 0) {
                Matcher mAttribute = patternAttribute.matcher(line);
                Matcher mData = patternData.matcher(line);

                //See if reached @DATA token
                if (mData.find()) {
                    dataFound = 1;

                } else if (mAttribute.find()) {
                    numAttribute++;

                    //If so, check if attribute is hierarchical, numeric or categoric
                    Matcher mNumeric = patternNumeric.matcher(line);
                    Matcher mHierarchical = patternHierarchical.matcher(line);

                    if (mHierarchical.find()) {
                        //Build structure to store the classes
                        setTreeClasses(line, "hierarchical");
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
    } catch (IOException ioe) {
        ioe.printStackTrace();
    }
}

public ArrayList<Integer> getPosClasses(String actualClasses) {
    ArrayList<Integer> positions = new ArrayList<Integer>();
    String[] vectorClasses = actualClasses.split("@");
    ArrayList<String> allClasses = new ArrayList<String>();

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

    HashSet hs = new HashSet();
    hs.addAll(allClasses);
    allClasses.clear();
    allClasses.addAll(hs);

    for (int i = 0; i < allClasses.size(); i++) {
        for (int j = 0; j < classes.length; j++) {
            if (allClasses.get(i).equals(classes[j]) == true) {
                positions.add(j);
                break;
            }
        }
    }

    return positions;
}

public void buildClassesStructureTrain() {
    int numberClasses = this.numberOfClasses;
    binaryClassesTrain = new ArrayList<int[]>();

    int[] binaryVector = new int[numberClasses];
    String actualClasses;
    ArrayList<Integer> posClasses;

    for (int i = 0; i < datasetTrain.size(); i++) {
        actualClasses = datasetTrain.get(i)[datasetTrain.get(i).length - 1];
        posClasses = getPosClasses(actualClasses);

        for (int j = 0; j < posClasses.size(); j++) {
            binaryVector[posClasses.get(j)] = 1;
        }

        binaryClassesTrain.add(binaryVector.clone());

        for (int j = 0; j < binaryVector.length; j++) {
            binaryVector[j] = 0;
        }
    }
}

public void buildClassesStructureTest() {
    int numberClasses = this.numberOfClasses;
    binaryClassesTest = new ArrayList<int[]>();

    int[] binaryVector = new int[numberClasses];
    String actualClasses;
    ArrayList<Integer> posClasses;

    for (int i = 0; i < datasetTest.size(); i++) {
        actualClasses = datasetTest.get(i)[datasetTest.get(i).length - 1];
        posClasses = getPosClasses(actualClasses);

        for (int j = 0; j < posClasses.size(); j++) {
            binaryVector[posClasses.get(j)] = 1;
        }

        binaryClassesTest.add(binaryVector.clone());

        for (int j = 0; j < binaryVector.length; j++) {
            binaryVector[j] = 0;
        }
    }
}

private void setTreeClassesWeights(double[] weightingScheme, String rootClass, ArrayList<Integer> superClassesPositionsAux) {
    rootClass = rootClass.concat("/[0-9]+");
    Pattern pattern = Pattern.compile(rootClass + "$");

    ArrayList<Integer> superClassesPositions = new ArrayList<Integer>();
    int numParents = 0;
    double sum = 0.0;
    Matcher m;

    for (int i = 0; i < classes.length; i++) {
        m = pattern.matcher(classes[i]);

        if (m.find()) {
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

private void setWeightingScheme(double[] weightingScheme) {
    ArrayList<Integer> topLevelClassesPositions = new ArrayList<Integer>();

    for (int i = 0; i < classes.length; i++) {
        if (classes[i].contains("/") == false) {
            topLevelClassesPositions.add(i);
        }
    }

    ArrayList<Integer> superClassesPositions = new ArrayList<Integer>();
    String rootClass = "";

    for (int i = 0; i < topLevelClassesPositions.size(); i++) {
        weightingScheme[topLevelClassesPositions.get(i)] = 0.75;
        rootClass = "^" + classes[topLevelClassesPositions.get(i)];
        superClassesPositions.add(topLevelClassesPositions.get(i));
        setTreeClassesWeights(weightingScheme, rootClass, superClassesPositions);
    }
}

public void setTreeClasses(String lineClasses, String tokenHierarchical) {
    String[] vetLine = lineClasses.split(tokenHierarchical);
    classes = vetLine[1].split(",");
    classes[0] = classes[0].trim();

    weightingScheme = new double[classes.length];
    setWeightingScheme(weightingScheme);

    ArrayList<ArrayList<Integer>> positionClassesLevels = new ArrayList<ArrayList<Integer>>();

    String rootClass = "";
    ArrayList<Integer> positions = new ArrayList<Integer>();
    Pattern pattern;
    Matcher m;

    for (int i = 0; i < Parameters.getNumLevels(); i++) {
        rootClass = rootClass.concat("[0-9]+");
        pattern = Pattern.compile("^" + rootClass + "$");

        for (int j = 0; j < classes.length; j++) {
            m = pattern.matcher(classes[j]);

            if (m.find()) {
                positions.add(j);
            }
        }

        rootClass = rootClass.concat("/");
        positionClassesLevels.add(positions);
    }

    //Set the position of the classes by level
    positionClassesLevel = positionClassesLevels;
}

public void setMeanClassLabelVectorAllTrain(String trainDatasetFile) {
    readTrainData(trainDatasetFile);
    buildClassesStructureTrain();
    meanClassLabelVectorAllClassesTrain = new double[classes.length];

    for (int i = 0; i < binaryClassesTrain.size(); i++) {
        for (int j = 0; j < classes.length; j++) {
            meanClassLabelVectorAllClassesTrain[j] += binaryClassesTrain.get(i)[j];
        }
    }

    for (int i = 0; i < meanClassLabelVectorAllClassesTrain.length; i++) {
        meanClassLabelVectorAllClassesTrain[i] = meanClassLabelVectorAllClassesTrain[i] / binaryClassesTrain.size();
    }
}

public void setMeanClassLabelVectorAllTest(String testDatasetFile) {
    readTestData(testDatasetFile);
    buildClassesStructureTest();
    meanClassLabelVectorAllClassesTest = new double[classes.length];
    for (int i = 0; i < binaryClassesTest.size(); i++) {
        for (int j = 0; j < classes.length; j++) {
            meanClassLabelVectorAllClassesTest[j] += binaryClassesTest.get(i)[j];
        }
    }

    for (int i = 0; i < meanClassLabelVectorAllClassesTest.length; i++) {
        meanClassLabelVectorAllClassesTest[i] = meanClassLabelVectorAllClassesTest[i] / binaryClassesTest.size();
    }
}

public ArrayList<String> getRuleAntecedents(String rule) {
    ArrayList<String> result = new ArrayList<>();
    int index;

    rule = rule.trim();
    index = rule.indexOf("=");
    rule = rule.substring(index + 1, rule.length());
    index = rule.indexOf("THEN");
    rule = rule.substring(0, index);

    if (!rule.contains("AND")) {
        result.add(rule);

    } else {
        while (!rule.isEmpty()) {
            index = rule.indexOf("AND");

            if (index == -1) {
                result.add(rule);
                break;
            }

            result.add(rule.substring(0, index));
            rule = rule.substring(index + 3, rule.length());
        }
    }

    return result;
}

public double[] getRuleConsequent(String rule) {
    int indexThen = rule.indexOf("THEN");
    String[] strSplit = rule.substring(indexThen + 5).split(" ");
    double[] result = new double[strSplit.length];

    for (int i = 0; i < strSplit.length; i++) {
        result[i] = Double.parseDouble(strSplit[i]);
    }

    return result;
}

public double[] getRuleConsequentClus(String rule) {
    rule = rule.replace("[", "");
    rule = rule.replace("]", "");

    int indexThen = rule.indexOf("THEN");
    String[] strSplit = rule.substring(indexThen + 5).split(",");
    double[] result = new double[strSplit.length];

    for (int i = 0; i < strSplit.length; i++) {
        result[i] = Double.parseDouble(strSplit[i]);
    }

    return result;
}

public int getAttrIndexTestedInARule(String antecedent) {
    antecedent = antecedent.trim();
    int index;
    String attr;

    if (!String.valueOf(antecedent.charAt(0)).matches("^[0-9]+$") && !antecedent.startsWith("-")) { // so is: attributeValue <= supLim  or attributeValue >= infLim, or > or <     
        if (antecedent.contains(">") && !antecedent.contains("=")) {
            index = antecedent.indexOf(">");
            attr = antecedent.substring(0, index);
            attr = attr.trim();

            return searchForAttr(attr);

        } else if (antecedent.contains("<") && !antecedent.contains("=")) {
            index = antecedent.indexOf("<");
            attr = antecedent.substring(0, index);
            attr = attr.trim();

            return searchForAttr(attr);

        } else { //<= or >=
            index = antecedent.indexOf("=");
            attr = antecedent.substring(0, index - 1);
            attr = attr.trim();

            return searchForAttr(attr);
        }

    } else { // then: // infLim <= attributeValue <= supLim
        index = antecedent.indexOf("=");
        attr = antecedent.substring(index + 2, antecedent.lastIndexOf("<="));
        attr = attr.trim();

        return searchForAttr(attr);
    }
}

public int searchForAttr(String attr) {
    for (int j = 0; j < attributes.length; j++) {
        if (attributes[j].contentEquals(attr)) {
            return j;
        }
    }

    return -1;
}

public int classifyWithARule(String[] example, String rule) {
    ArrayList<String> antecedents = getRuleAntecedents(rule);
//    ArrayList<String> antecedents = getRuleAntecedentsClus(rule);
    int coverage = 0;
    int index;
    Double infLim;
    Double supLim;
    int posAttribute;
    double attributeValue;

    for (int i = 0; i < antecedents.size(); i++) {
        String test = antecedents.get(i);
        test = test.trim();
        posAttribute = getAttrIndexTestedInARule(test);
        attributeValue = Double.parseDouble(example[posAttribute]);

        if (test.contains(">=")) {
            index = test.indexOf(">=");
            infLim = Double.parseDouble(test.substring(index + 2, test.length()));
            coverage = greaterEqual(attributeValue, infLim);

        } else if (test.startsWith("-") || String.valueOf(test.charAt(0)).matches("^[0-9]+$")) { //<= att <=
            infLim = Double.parseDouble(test.substring(0, test.indexOf("<=")));
            supLim = Double.parseDouble(test.substring(test.lastIndexOf("<=") + 2, test.length()));
            coverage = compoundTerm(infLim, supLim, attributeValue);

        } else if (test.contains("<=")) {
            index = test.indexOf("<=");
            supLim = Double.parseDouble(test.substring(index + 2, test.length()));
            coverage = lessEqual(attributeValue, supLim);

        } else if (test.contains(">")) {
            index = test.indexOf(">");
            infLim = Double.parseDouble(test.substring(index + 1, test.length()));
            coverage = greater(attributeValue, infLim);

        } else {
            index = test.indexOf("<");
            supLim = Double.parseDouble(test.substring(index + 1, test.length()));
            coverage = less(attributeValue, supLim);

        }

        if (coverage == 0) {
            return coverage;
        }
    }

    return coverage;
}

private static int lessEqual(double attributeValue, double supLim) {
    if (attributeValue <= supLim) {
        return 1;
    } else {
        return 0;
    }
}

private static int less(double attributeValue, double supLim) {
    if (attributeValue < supLim) {
        return 1;
    } else {
        return 0;
    }
}

public static int greaterEqual(double attributeValue, double infLim) {
    if (attributeValue >= infLim) {
        return 1;
    } else {
        return 0;
    }
}

public static int greater(double attributeValue, double infLim) {
    if (attributeValue > infLim) {
        return 1;
    } else {
        return 0;
    }
}

private static int compoundTerm(double infLim, double supLim, double attributeValue) {
    int left = greaterEqual(attributeValue, infLim);
    int right = lessEqual(attributeValue, supLim);

    if (left == 1 && right == 1) {
        return 1;
    } else {
        return 0;
    }
}

public String[] leConfig() {
    String regExp[] = {"nameDataset =",
        "specificFold =",
        "folds =",
        "pathTrainDataset =",
        "pathTestDataset =",
        "pathRules =",
        "pathToSavePredictions ="};

    Pattern comment = Pattern.compile("#");
    String[] result = new String[7];

    for (int i = 0; i < regExp.length; i++) {
        try {
//                FileReader reader = new FileReader("/home/geantrindade/Dropbox/posGrad/GACerriMaven/src/main/java/hmc_ga/teste.txt");
//            FileReader reader = new FileReader("teste.txt");
            FileReader reader = new FileReader("C:\\Users\\gean_\\Dropbox\\posGrad\\GAs\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\config\\test.txt");
            BufferedReader buffReader = new BufferedReader(reader);

            Pattern pattern = Pattern.compile(regExp[i]);
            String line = null;

            while ((line = buffReader.readLine()) != null) {
                Matcher m = pattern.matcher(line);
                Matcher m1 = comment.matcher(line);

                if (m.find() && !m1.find()) {
                    String[] vectorLine = line.split(" = ");
                    result[i] = vectorLine[1];
                    break;
                }
            }

            buffReader.close();
            reader.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }
    return result;
}

public void makePredictions(String rulesFile, String fullPathToSavePrediction) {
    ArrayList<String> rules = readRulesFile(rulesFile);
//    ArrayList<String> rules = readRulesFileClus(rulesFile);
    String[] example;
    int coverage;
    double[][] matrixPredictions = new double[binaryClassesTest.size()][binaryClassesTest.get(0).length];

    for (int i = 0; i < datasetTest.size(); i++) {
        example = datasetTest.get(i);
        coverage = 0;

        for (int j = 0; j < rules.size(); j++) {
            String rule = rules.get(j);
            coverage = classifyWithARule(example, rule);

            if (coverage == 1) {
                double[] vectorConsequent = getRuleConsequent(rule);
//                double[] vectorConsequent = getRuleConsequentClus(rule);
                System.arraycopy(vectorConsequent, 0, matrixPredictions[i], 0, vectorConsequent.length);
                break;
            }
        }

        //If no rule classify the example, apply the default rule  
        if (coverage == 0) {
            double[] defaultRule = meanClassLabelVectorAllClassesTrain.clone();
            System.arraycopy(defaultRule, 0, matrixPredictions[i], 0, defaultRule.length);
        }
    }

    try {
        PrintWriter writer = new PrintWriter(fullPathToSavePrediction + "originalPredictions.txt", "UTF-8");

        for (int i = 0; i < matrixPredictions.length; i++) {
            for (int j = 0; j < matrixPredictions[i].length; j++) {
                writer.print(matrixPredictions[i][j] + " ");
            }

            writer.print("\n");
        }

        writer.flush();
        writer.close();

    } catch (FileNotFoundException ex) {
        Logger.getLogger(Validation.class
                .getName()).log(Level.SEVERE, null, ex);

    } catch (IOException ex) {
        Logger.getLogger(ClassificationWithRule.class.getName()).log(Level.SEVERE, null, ex);
    }
}

public double fMeasurePrediction(int indexExample, double[] prediction) {
    double fmeasure = 0;
    double sumIntersection = 0;
    double minSumPredicted = 0;
    double sumReal = 0;
    double sumPredicted = 0;

    for (int j = 0; j < classes.length; j++) {
//        if (prediction[j] >= 0.5 && binaryClassesTest.get(indexExample)[j] == 1) {
        if (prediction[j] >= 0.5 && binaryClassesTrain.get(indexExample)[j] == 1) {
            sumIntersection++;
        }

        if (prediction[j] > 0) {
            sumPredicted++;
        }

//        if (binaryClassesTest.get(indexExample)[j] == 1) {
        if (binaryClassesTrain.get(indexExample)[j] == 1) {
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

//for all rules
public ArrayList<ArrayList<Integer>> obtainRulesIndexCoveredExamples(ArrayList<String> rules) {
    ArrayList<ArrayList<Integer>> indexCoveredExamples = new ArrayList<>();

    for (int i = 0; i < rules.size(); i++) {
        indexCoveredExamples.add(i, new ArrayList<Integer>());
    }

    for (int i = 0; i < datasetTrain.size(); i++) {
        String[] example = datasetTrain.get(i);

        for (int j = 0; j < rules.size(); j++) {
            int coverage = classifyWithARule(example, rules.get(j));

            if (coverage == 1) { //Rule satisfies example
                indexCoveredExamples.get(j).add(i);
            }
        }
    }

    return indexCoveredExamples;
}

//just rules who can cover the test example
public int getExamplesBasedOnSimilarity(String[] testExample, ArrayList<String> selectedRules, ArrayList<ArrayList<Integer>> selectedRulesIndexCoveredExamples, int rangeOfConsideration) {
    int result = 0;

    for (int i = 0; i < selectedRules.size(); i++) {
        if (checkSimilarity(testExample, selectedRulesIndexCoveredExamples.get(i), rangeOfConsideration)) {
            result += 1;
        }
    }

    return result;
}

//return the indexes
public ArrayList<Integer> getExamplesBasedOnSimilarity(String[] testExample, int rangeOfConsideration) {
    ArrayList<Integer> trainExamples = new ArrayList<>();

    for (int i = 0; i < datasetTrain.size(); i++) {
        if (checkSimilarityBetweenExamples(testExample, datasetTrain.get(i), rangeOfConsideration)) {
            trainExamples.add(i);
        }
    }

    return trainExamples;
}

public Integer getMostSimilarExampleBasedOnSum(ArrayList<Integer> similarTrainExamples, String[] testExample) {
    int valueTest = 0;
    for (int i = 0; i < testExample.length - 1; i++) {
        valueTest += Integer.valueOf(testExample[i]);
    }

    int valueTrain = 0, minDiferrence = 0, diferrenceAux = 0, indexMostSimilar = 0;
    for (int i = 0; i < similarTrainExamples.size(); i++) {
        String[] trainExample = datasetTrain.get(similarTrainExamples.get(i));

        for (int j = 0; j < trainExample.length - 1; j++) {
            valueTrain += Integer.valueOf(trainExample[j]);
        }

        diferrenceAux = Math.abs(valueTest - valueTrain);

        if (minDiferrence == 0) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        if (diferrenceAux < minDiferrence) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        valueTrain = 0;
        diferrenceAux = 0;
    }

    return indexMostSimilar;
}

public Integer getMostSimilarExampleBasedOnAverage(ArrayList<Integer> similarTrainExamples, String[] testExample) {
    int valueTest = 0;
    for (int i = 0; i < testExample.length - 1; i++) {
        valueTest += Integer.valueOf(testExample[i]);
    }
    valueTest = valueTest / testExample.length;

    int valueTrain = 0, minDiferrence = 0, diferrenceAux = 0, indexMostSimilar = 0;
    for (int i = 0; i < similarTrainExamples.size(); i++) {
        String[] trainExample = datasetTrain.get(similarTrainExamples.get(i));

        for (int j = 0; j < trainExample.length - 1; j++) {
            valueTrain += Integer.valueOf(trainExample[j]);
        }
        valueTrain = valueTrain / trainExample.length;

        diferrenceAux = Math.abs(valueTest - valueTrain);

        if (minDiferrence == 0) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        if (diferrenceAux < minDiferrence) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        valueTrain = 0;
        diferrenceAux = 0;
    }

    return indexMostSimilar;
}

public Integer getMostSimilarExampleBasedOnEachDifference(ArrayList<Integer> similarTrainExamples, String[] testExample) {
    int[] valuesTest = new int[testExample.length];

    for (int i = 0; i < testExample.length - 1; i++) {
        valuesTest[i] = Integer.valueOf(testExample[i]);
    }

    int[] valuesTrain = new int[testExample.length];
    int minDiferrence = 0, diferrenceAux = 0, indexMostSimilar = 0;

    for (int i = 0; i < similarTrainExamples.size(); i++) {
        String[] trainExample = datasetTrain.get(similarTrainExamples.get(i));

        for (int j = 0; j < trainExample.length - 1; j++) {
            valuesTrain[j] = Integer.valueOf(trainExample[j]);
            diferrenceAux += Math.abs(valuesTest[j] - valuesTrain[j]); //removes minus signal
        }

        if (minDiferrence == 0) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        if (diferrenceAux < minDiferrence) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        diferrenceAux = 0;
    }

    return indexMostSimilar;
}

public Integer getMostSimilarExampleBasedOnEachDifferenceAverage(ArrayList<Integer> similarTrainExamples, String[] testExample) {
    int[] valuesTest = new int[testExample.length];

    for (int i = 0; i < testExample.length - 1; i++) {
        valuesTest[i] = Integer.valueOf(testExample[i]);
    }

    int[] valuesTrain = new int[testExample.length];
    int minDiferrence = 0, diferrenceAux = 0, indexMostSimilar = 0;

    for (int i = 0; i < similarTrainExamples.size(); i++) {
        String[] trainExample = datasetTrain.get(similarTrainExamples.get(i));

        for (int j = 0; j < trainExample.length - 1; j++) {
            valuesTrain[j] = Integer.valueOf(trainExample[j]);
            diferrenceAux += Math.abs(valuesTest[j] - valuesTrain[j]); //removes minus signal
        }

        diferrenceAux = diferrenceAux / trainExample.length - 1;

        if (minDiferrence == 0) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        if (diferrenceAux < minDiferrence) {
            minDiferrence = diferrenceAux;
            indexMostSimilar = i;
        }

        diferrenceAux = 0;
    }

    return indexMostSimilar;
}

// checar quantas instancias tem similaridade, ai pegar a regra com maior numero
public boolean checkSimilarityBetweenExamples(String[] testExample, String[] trainExample, int rangeOfConsideration) {
    // until 40% of similarity, lower than that it is not acceptable (134/336) = 202
    // until 50% of similarity, lower than that it is not acceptable (168/336) = 168
    // until 60% of similarity, lower than that it is not acceptable (202/336) = 134
    // until 70% of similarity, lower than that it is not acceptable (235/336) = 101
    // until 80% of similarity, lower than that it is not acceptable (269/336) = 67
    // until 90% of similarity, lower than that it is not acceptable (302/336) = 34
    // until 95% of similarity, lower than that it is not acceptable (319/336) = 17
    // until 98% of similarity, lower than that it is not acceptable (329/336) = 7
    int tolerance = 0;

    for (int j = 0; j < trainExample.length - 1; j++) {
        int valueTrain = Integer.valueOf(trainExample[j]);
        int valueTest = Integer.valueOf(testExample[j]);

        if (valueTest >= (valueTrain - rangeOfConsideration) && valueTest <= (valueTrain + rangeOfConsideration)) {
            if (j == trainExample.length - 2) {
                return true;
            }

        } else if (tolerance < 34) {
            if (j == trainExample.length - 2) {
                return true;
            }
            tolerance++;

        } else {
            j = trainExample.length;
        }
    }

    return false;
}

// checar quantas instancias tem similaridade, ai pegar a regra com maior numero
public boolean checkSimilarity(String[] testExample, ArrayList<Integer> indexCoveredExamples, int rangeOfConsideration) {
    for (int i = 0; i < indexCoveredExamples.size(); i++) {
        String[] trainExample = datasetTrain.get(indexCoveredExamples.get(i));

        for (int j = 0; j < trainExample.length - 1; j++) {
            int valueTrain = Integer.valueOf(trainExample[j]);
            int valueTest = Integer.valueOf(testExample[j]);

            if (valueTest >= valueTrain - rangeOfConsideration && valueTest <= valueTrain + rangeOfConsideration) {
                if (j == trainExample.length - 2) {
                    return true;
                }
            } else {
                j = trainExample.length;
            }
        }
    }

    return false;
}

public int getBestRuleBasedOnCountMax(ArrayList<String> selectedRules) {
    for (int i = 0; i < selectedRules.size(); i++) {
        double[] consequent = getRuleConsequent(selectedRules.get(i));

        if (checkMaxs(consequent)) {
            return i;
        }
    }

    return 0;
}

// checar quantas instancias tem similaridade, ai pegar a regra com maior numero
public boolean checkMaxs(double[] ruleConsequent) {
    return ruleConsequent[0] == 1.0 || ruleConsequent[6] == 1.0;
}

public int indexBestRuleFromAGroup(ArrayList<String> rules) {
    ArrayList<String> rulesGroup1 = new ArrayList<>();
    ArrayList<String> rulesGroup2 = new ArrayList<>();

    for (int i = 0; i < rules.size(); i++) {
        double[] consequent = getRuleConsequent(rules.get(i));

        if (consequent[0] == 1.0) {
            rulesGroup1.add(rules.get(i));
        } else {
            rulesGroup2.add(rules.get(i));
        }
    }

    return 0;
}

public String getBestRuleBasedOnGroupAndSimilarTrainExamples(ArrayList<Integer> trainExamplesIndexes,
        ArrayList<String> rulesGroup1, ArrayList<String> rulesGroup2) {

    if (trainExamplesIndexes.isEmpty() || rulesGroup1.isEmpty() || rulesGroup2.isEmpty()) {
        return rulesGroup1.get(0);

    } else {
        double fMeanG1 = 0.0, fMeanG2 = 0.0, temp = 0.0;

        for (int i = 0; i < trainExamplesIndexes.size(); i++) {
            for (int j = 0; j < rulesGroup1.size(); j++) {
                double[] consequent = getRuleConsequent(rulesGroup1.get(j));
                temp += fMeasurePrediction(trainExamplesIndexes.get(i), consequent);
            }
            fMeanG1 += temp / rulesGroup1.size();
            temp = 0.0;

            for (int j = 0; j < rulesGroup2.size(); j++) {
                double[] consequent = getRuleConsequent(rulesGroup2.get(j));
                temp += fMeasurePrediction(trainExamplesIndexes.get(i), consequent);
            }
            fMeanG2 += temp / rulesGroup2.size();
            temp = 0.0;
        }

        if (fMeanG1 >= fMeanG2) {
            return rulesGroup1.get(0); //first rule of the set 1...
        } else {
            return rulesGroup2.get(0); //firs rule of the set 2...
        }
    }
}

public void makePredictionsImproved(String pathRules, String fullPathToSavePrediction) {
    ArrayList<String> rules = readRulesFile(pathRules);
    ArrayList<ArrayList<String>> examplesAndRules = new ArrayList<>();
    String[] example;
    HashSet<String> rulesOnceUsed = new HashSet<>();
    int coverage;

    double[][] matrixPredictions = new double[binaryClassesTest.size()][binaryClassesTest.get(0).length];

    for (int i = 0; i < datasetTest.size(); i++) {
        example = datasetTest.get(i);
        examplesAndRules.add(i, new ArrayList<String>());

        for (int j = 0; j < rules.size(); j++) {
            String rule = rules.get(j);
            coverage = classifyWithARule(example, rule);

            if (coverage == 1) {
                examplesAndRules.get(i).add(rule);
            }
        }

        //passar thr nos consequentes das regras pra pegar 1s, contar 1s, e dar pesos pras regras, da primeira a ultima... a que tirar maior pontuacao vence!
        if (!examplesAndRules.get(i).isEmpty()) {
            //store the indexes
            ArrayList<Integer> trainExamplesSimilar = new ArrayList<>();
            int range = 0;

//            while (trainExamplesSimilar.isEmpty() && range < 201) {
            while (trainExamplesSimilar.isEmpty()) {
                trainExamplesSimilar = getExamplesBasedOnSimilarity(example, range);
                range++;
            }
//            System.out.println("range: " + (range - 1));

//            double temp, betterValue = 0.0;
//            int betterRuleIndex = 0;
//
//            for (int j = 0; j < examplesAndRules.get(i).size(); j++) {
//                double[] consequent = getRuleConsequent(examplesAndRules.get(i).get(j));
//                    temp = getConsequentSum(consequent);
//                    temp = getConsequentAverage(consequent);
//                    temp = getNumberOfActivatedClasses(consequent);
//                    temp = getActivatedClassesMultiLevels(consequent);
//
//                if (trainExamplesSimilar.isEmpty()) {
//                    j = examplesAndRules.get(i).size();
//                    j = 0;
//
//                } else {
//                    for (int k = 0; k < trainExamplesSimilar.size(); k++) {
//                        int exampleIndex = trainExamplesSimilar.get(k);
//                        temp += fMeasurePrediction(i, consequent);
//                        temp = fMeasurePrediction(i, consequent);
//
//                    }
//                    if ((temp / trainExamplesSimilar.size()) > betterValue) {
//                    if (temp > betterValue) {
//                        betterValue = temp;
//                        betterRuleIndex = j;
//                    } else if (temp == betterValue) {
//                        double[] bestConsequent = getRuleConsequent(examplesAndRules.get(i).get(betterRuleIndex));
//                        double bestTemp = getNumberOfActivatedClasses(bestConsequent);
//
//                        temp = getConsequentSum(consequent);
//
//                        if (temp >= bestTemp) {
//                            betterRuleIndex = j;
//                        }
//                    }
//
//                    temp = 0.0;
//                }
//            }
            int betterRuleIndex = 0;
//            boolean defaultRuleWasUsed = false;

            if (!trainExamplesSimilar.isEmpty()) {
//                HashMap<Integer, String> retrotransposonRules = new HashMap<>();
//                HashMap<Integer, String> transposonRules = new HashMap<>();

                double temp = 0.0, betterValue = 0.0;

//                System.out.println("number of rules: "+examplesAndRules.get(i).size());
                
//                System.out.println("trainExamplesSimilar size: " + trainExamplesSimilar.size());

                for (int j = 0; j < examplesAndRules.get(i).size(); j++) {
                    String rule = examplesAndRules.get(i).get(j);
                    double[] consequent = getRuleConsequent(rule);

//                    for (int k = 0; k < consequent.length; k++) {
//                        System.out.print(consequent[k] + " ");
//                    }
//                    System.out.println("");
//                    if (consequent[0] >= 0.5) { //equals to "1" class
//                        retrotransposonRules.put(j, rule);
//                    } else if (consequent[6] >= 0.5) {//equals to "2" class
//                        transposonRules.put(j, rule);
//                    }
//                }
//                double fMeasureRetroAverage = 0.0;
//                    Random rn = new Random();
//                    int index = rn.nextInt(trainExamplesSimilar.size());
//int index = getMostSimilarExampleBasedOnEachDifference(trainExamplesSimilar, example);
//                    temp = fMeasurePrediction(trainExamplesSimilar.get(index), consequent);
                    for (int k = 0; k < trainExamplesSimilar.size(); k++) {
                        temp += fMeasurePrediction(trainExamplesSimilar.get(k), consequent);
                    }

//                    System.out.println("temp " + j + ": " + temp);

                    temp = temp / trainExamplesSimilar.size();

//                    System.out.println("temp average " + j + ": " + temp);
                    
                    if (temp > betterValue) {
                        betterValue = temp;
                        betterRuleIndex = j;
                    }

                    temp = 0.0;
//                double tempR, betterValueR = 0.0;
//                Integer keyBestRulesGroupR = 0;
//                for (Integer key : retrotransposonRules.keySet()) {
//                    double[] consequent = getRuleConsequent(retrotransposonRules.get(key));
//
//                    tempR = fMeasurePrediction(trainExamplesSimilar.get(index), consequent);
//                    fMeasureRetroAverage += tempR;
//
//                    if (tempR > betterValueR) {
//                        betterValueR = tempR;
//                        keyBestRulesGroupR = key;
//                    }
//                }
//                if (!retrotransposonRules.isEmpty() || fMeasureRetroAverage == 0.0) {
//                    fMeasureRetroAverage = fMeasureRetroAverage / retrotransposonRules.size();
//                }
//                double fMeasureTransAverage = 0.0;
//                double tempT, betterValueT = 0.0;
//                Integer keyBestRulesGroupT = 0;
//
//                for (Integer key : transposonRules.keySet()) {
//                    double[] consequent = getRuleConsequent(transposonRules.get(key));
//
//                    tempT = fMeasurePrediction(trainExamplesSimilar.get(index), consequent);
//                    fMeasureTransAverage += tempT;
//
//                    if (tempT > betterValueT) {
//                        betterValueT = tempT;
//                        keyBestRulesGroupT = key;
//                    }
//                }
//
//                if (!transposonRules.isEmpty() || fMeasureTransAverage == 0.0) {
//                    fMeasureTransAverage = fMeasureTransAverage / transposonRules.size();
//                }
                    //se nao houver instancias semelhantes com range menor que 5, ver se existem mais regras trans ou retro e escolher a primeira do conjunto
//                if (fMeasureRetroAverage >= fMeasureTransAverage) {
//                    betterRuleIndex = keyBestRulesGroupR;
//
//                } else if (fMeasureRetroAverage < fMeasureTransAverage) {
//                    betterRuleIndex = keyBestRulesGroupT;
//                }
                }

//                if (betterValue == 0.0) {
//                    System.out.println("default rule cause bad fmeasure happened");
//                    //If no rule classify the example, apply the default rule
//                    double[] defaultRule = meanClassLabelVectorAllClassesTrain.clone();
//                    System.arraycopy(defaultRule, 0, matrixPredictions[i], 0, defaultRule.length);
//                    defaultRuleWasUsed = true;
//                    Random rn = new Random();
//                    int index = rn.nextInt(examplesAndRules.get(i).size());
//                    betterRuleIndex = index;
//                }
            } //else if (trainExamplesSimilar.isEmpty() && !rulesOnceUsed.isEmpty()) {
//                double temp, betterValue = 0.0;
//                
//                for (int j = 0; j < examplesAndRules.get(i).size(); j++) {
//                    if (rulesOnceUsed.contains(examplesAndRules.get(i).get(j))) {
////                        betterRuleIndex = j;
////                        j = examplesAndRules.get(i).size();
//                    
//double[] consequent = getRuleConsequent(examplesAndRules.get(i).get(j));
//                    temp = fMeasurePrediction(example, consequent);
//                    fMeasureTransAverage += tempT;
//
//
//                        System.out.println("regra ja usada!");
//                    }
//                }
//            }
            System.out.println("betterRuleIndex: " + betterRuleIndex);
//            if (defaultRuleWasUsed == false) {
            double[] vectorConsequent = getRuleConsequent(examplesAndRules.get(i).get(betterRuleIndex));
            System.arraycopy(vectorConsequent, 0, matrixPredictions[i], 0, vectorConsequent.length);
//            }
//            rulesOnceUsed.add(examplesAndRules.get(i).get(betterRuleIndex));

//            System.out.print("consq: ");
//            for (int j = 0; j < vectorConsequent.length; j++) {
//                System.out.print(vectorConsequent[j] + " ");
//            }
//            System.out.println("\n");
//                double[] vectorConsequent = getRuleConsequentClus(rulesAndExamples.get(i).get(betterRuleIndex));
        } else {
            System.out.println("default rule");
            //If no rule classify the example, apply the default rule
            double[] defaultRule = meanClassLabelVectorAllClassesTrain.clone();
            System.arraycopy(defaultRule, 0, matrixPredictions[i], 0, defaultRule.length);
        }
    }

//    System.out.println("rules once used size: " + rulesOnceUsed.size());
    try {
        PrintWriter writer = new PrintWriter(fullPathToSavePrediction + "improvedPredictions.txt", "UTF-8");

        for (int i = 0; i < matrixPredictions.length; i++) {
            for (int j = 0; j < matrixPredictions[i].length; j++) {
                writer.print(matrixPredictions[i][j] + " ");
            }

            writer.print("\n");
        }

        writer.flush();
        writer.close();

    } catch (FileNotFoundException ex) {
        Logger.getLogger(Validation.class
                .getName()).log(Level.SEVERE, null, ex);

    } catch (UnsupportedEncodingException ex) {
        Logger.getLogger(Validation.class
                .getName()).log(Level.SEVERE, null, ex);
    }
}

public HashMap getRankingOfRulesByFmeasure(ArrayList<String> rules) {
    HashMap<String, Double> fitness = new HashMap<>();
    ArrayList<Integer> numberCoveredExamples = new ArrayList<>();
    ArrayList<ArrayList<Integer>> indexCoveredExamples = new ArrayList<>();
    ArrayList<double[]> meanClassLabelVectorCovered = new ArrayList<>();

    for (int i = 0; i < rules.size(); i++) {
        numberCoveredExamples.add(0);
        indexCoveredExamples.add(i, new ArrayList<Integer>());
        meanClassLabelVectorCovered.add(new double[classes.length]);
    }

    for (int i = 0; i < datasetTrain.size(); i++) {
        String[] example = datasetTrain.get(i);

        for (int j = 0; j < rules.size(); j++) {
            int coverage = classifyWithARule(example, rules.get(j));
            numberCoveredExamples.set(j, numberCoveredExamples.get(j) + coverage);

            if (coverage == 1) { //Rule satisfies example
                indexCoveredExamples.get(j).add(i);
            }
        }
    }

    ArrayList<int[]> binaryClasses = binaryClassesTrain;

    for (int i = 0; i < rules.size(); i++) {

        for (int j = 0; j < indexCoveredExamples.get(i).size(); j++) {
            int[] binaryVector = binaryClasses.get(indexCoveredExamples.get(i).get(j));

            for (int k = 0; k < binaryVector.length; k++) {
                meanClassLabelVectorCovered.get(i)[k] += binaryVector[k];
            }
        }

        if (indexCoveredExamples.get(i).size() > 0) {
            for (int j = 0; j < meanClassLabelVectorCovered.get(i).length; j++) {
                meanClassLabelVectorCovered.get(i)[j] = meanClassLabelVectorCovered.get(i)[j] / indexCoveredExamples.get(i).size();
            }
        }

        double[][] matrixPredictions = new double[numberCoveredExamples.get(i)][classes.length];
        for (int j = 0; j < numberCoveredExamples.get(i); j++) {
            System.arraycopy(meanClassLabelVectorCovered.get(i), 0, matrixPredictions[j], 0, classes.length);
        }

        fitness.put(rules.get(i), evaluationFmeasureFitness(matrixPredictions, indexCoveredExamples.get(i)));
    }

    return fitness;
}

public double evaluationFmeasureFitness(double[][] predictedClasses, ArrayList<Integer> indexExamples) {
    double fmeasure = 0;
    double sumIntersection = 0;
    double minSumPredicted = 0;
    double sumReal = 0;

    //Matrix to store the outputs on the test data
    int[][] binaryMatrix = new int[indexExamples.size()][classes.length];

    applyThresholds(binaryMatrix, predictedClasses, 0.5, 0);

    for (int i = 0; i < indexExamples.size(); i++) {
        int numInst = indexExamples.get(i);

        double sumPredictedExample = 0;
        double sumRealExample = 0;

        for (int j = 0; j < classes.length; j++) {
            if (binaryMatrix[i][j] == 1 && binaryClassesTrain.get(numInst)[j] == 1) {
                sumIntersection++;
            }

            if (binaryMatrix[i][j] == 1) {
                sumPredictedExample++;
            }

            if (binaryClassesTrain.get(numInst)[j] == 1) {
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

public HashMap getRankingOfRulesByFitness(ArrayList<String> rules) {
    HashMap<String, Double> fitness = new HashMap<>();
    ArrayList<Integer> numberCoveredExamples = new ArrayList<>();
    ArrayList<ArrayList<Integer>> indexCoveredExamples = new ArrayList<>();
    ArrayList<ArrayList<Integer>> indexUncoveredExamples = new ArrayList<>();
    ArrayList<double[]> meanClassLabelVectorCovered = new ArrayList<>();
    ArrayList<double[]> meanClassLabelVectorUncovered = new ArrayList<>();

    for (int i = 0; i < rules.size(); i++) {
        numberCoveredExamples.add(0);
        indexCoveredExamples.add(i, new ArrayList<Integer>());
        indexUncoveredExamples.add(i, new ArrayList<Integer>());
        meanClassLabelVectorCovered.add(new double[classes.length]);
        meanClassLabelVectorUncovered.add(new double[classes.length]);
    }

    for (int i = 0; i < datasetTrain.size(); i++) {
        String[] example = datasetTrain.get(i);

        for (int j = 0; j < rules.size(); j++) {
            int coverage = classifyWithARule(example, rules.get(j));
            numberCoveredExamples.set(j, numberCoveredExamples.get(j) + coverage);

            if (coverage == 1) { //Rule satisfies example
                indexCoveredExamples.get(j).add(i);
            } else {
                indexUncoveredExamples.get(j).add(i);
            }
        }
    }

    ArrayList<int[]> binaryClasses = binaryClassesTrain;

    for (int i = 0; i < rules.size(); i++) {
        //Covered examples
        for (int j = 0; j < indexCoveredExamples.get(i).size(); j++) {
            int[] binaryVector = binaryClasses.get(indexCoveredExamples.get(i).get(j));

            for (int k = 0; k < binaryVector.length; k++) {
                meanClassLabelVectorCovered.get(i)[k] += binaryVector[k];
            }
        }

        if (indexCoveredExamples.get(i).size() > 0) {
            for (int j = 0; j < meanClassLabelVectorCovered.get(i).length; j++) {
                meanClassLabelVectorCovered.get(i)[j] = meanClassLabelVectorCovered.get(i)[j] / indexCoveredExamples.get(i).size();
            }
        }

        //Uncovered examples
        for (int j = 0; j < indexUncoveredExamples.get(i).size(); j++) {
            int[] binaryVector = binaryClasses.get(indexUncoveredExamples.get(i).get(j));

            for (int k = 0; k < binaryVector.length; k++) {
                meanClassLabelVectorUncovered.get(i)[k] += binaryVector[k];
            }
        }

        if (indexUncoveredExamples.get(i).size() > 0) {
            for (int j = 0; j < meanClassLabelVectorUncovered.get(i).length; j++) {
                meanClassLabelVectorUncovered.get(i)[j] = meanClassLabelVectorUncovered.get(i)[j] / indexUncoveredExamples.get(i).size();
            }
        }
    }

    setMeanClassLabelVectorAllTrain(pathTrainDataset + nameDataset.toLowerCase() + specificFold + "trainatt.arff");

    ArrayList<Double> AUPRCs = getAUPRC(numberCoveredExamples, meanClassLabelVectorCovered, indexCoveredExamples);
    ArrayList<Double> VarianceGains = getVarianceGain(indexCoveredExamples, indexUncoveredExamples, numberCoveredExamples,
            meanClassLabelVectorCovered, meanClassLabelVectorUncovered);

    if (AUPRCs.size() == VarianceGains.size()) {
        for (int i = 0; i < rules.size(); i++) {
            fitness.put(rules.get(i), (0.3 * AUPRCs.get(i)) + (0.7 * VarianceGains.get(i)));
        }
    } else {
        System.out.println("fatal error");
        System.exit(0);
    }

    return fitness;
}

public void countRules(String rule) {
    int cont = -1, index = 0;

    while (index != -1) {
        cont++;
        index = rule.indexOf("THEN", index);
        rule = rule.replaceFirst("THEN", "");
//        if(index != -1){
//            rule = rule.substring(index+4, rule.length());
//            cont++;
//        }
    }
    numberOfRules = cont;
}

public double getConsequentSum(double[] rule) {
    double sum = 0;
    for (int i = 0; i < rule.length; i++) {
        sum += rule[i];
    }

    return sum;
}

public double getConsequentAverage(double[] rule) {
    return (getConsequentSum(rule) / rule.length);
}

public double getNumberOfActivatedClasses(double[] rule) {
    double cont = 0;

    for (int i = 0; i < rule.length; i++) {
        if (rule[i] >= 0.5) {
            cont++;
        }
    }

    return cont;
}

public double getActivatedClassesMultiLevels(double[] rule) {
    double value = 0;
    double temp = 0;
//    1, 1.1, 1.1.1, 1.1.2, 1.4, 1.5, 2, 2.1, 2.1.1, 2.1.1.1, 2.1.1.2, 2.1.1.3, 2.1.1.8, 2.1.1.9

    double[] weightingAux = {0.95, 0.75, 0.6, 0.6, 0.75, 0.75, 0.95, 0.75, 0.45, 0.3, 0.3, 0.3, 0.3, 0.3};
    for (int i = 0; i < rule.length; i++) {
        if (rule[i] >= 0.5) {
            temp = rule[i] * weightingAux[i];
            value += temp;
        }
    }

    return value;
}

public double getActivatedClassesMultiLevelsAndAverage(double[] rule) {
    double value = 0;
    double temp = 0;
    for (int i = 0; i < rule.length; i++) {
//        if (rule[i] >= 0.5) {
        temp = rule[i] * weightingScheme[i];
        value += temp;
//        }
    }

    return value / rule.length;
}

public void escreveArquivo(String pathWithFile, String subject) {
    FileWriter writer;
    BufferedWriter bufferWriter;

    try {
        writer = new FileWriter(pathWithFile);
        bufferWriter = new BufferedWriter(writer);

        writer.write(subject);
        writer.flush();
        writer.close();

    } catch (IOException ex) {
        Logger.getLogger(Validation.class
                .getName()).log(Level.SEVERE, null, ex);
    }
}

public void run(String prediction) {
    if (specificFold > 0) {
        setMeanClassLabelVectorAllTrain(pathTrainDataset + nameDataset.toLowerCase() + specificFold + "trainatt.arff");

        setMeanClassLabelVectorAllTest(pathTestDataset + nameDataset.toLowerCase() + specificFold + "testatt.arff");

        if (prediction.equals("makePredictions")) {
            makePredictions(pathRules + "rules.txt",
                    //            makePredictions(pathRules + "fold" + specificFold + "\\rules.txt",
                    //                    pathToSavePredictions + "fold" + specificFold + "\\");
                    pathToSavePredictions);

        } else if (prediction.equals("makePredictionsImproved")) {
            makePredictionsImproved(pathRules + "rules.txt",
                    //            makePredictionsImproved(pathRules + "fold" + specificFold + "\\rules.txt",
                    //                    pathToSavePredictions + "fold" + specificFold + "\\");
                    pathToSavePredictions);
        }

    } else if (specificFold == -1) {
        if (prediction.equals("makePredictions")) {
            for (int i = 1; i <= folds; i++) {
                setMeanClassLabelVectorAllTrain(pathTrainDataset + nameDataset.toLowerCase() + i + "trainatt.arff");
                setMeanClassLabelVectorAllTest(pathTestDataset + nameDataset.toLowerCase() + i + "testatt.arff");
                makePredictions(pathRules + "fold" + i + "\\rules.txt",
                        pathToSavePredictions + "fold" + i + "\\");
            }

        } else if (prediction.equals("makePredictionsImproved")) {
            for (int i = 1; i <= folds; i++) {
                setMeanClassLabelVectorAllTrain(pathTrainDataset + nameDataset.toLowerCase() + i + "trainatt.arff");
                setMeanClassLabelVectorAllTest(pathTestDataset + nameDataset.toLowerCase() + i + "testatt.arff");
                makePredictionsImproved(pathRules + "fold" + i + "\\rules.txt",
                        pathToSavePredictions + "fold" + i + "\\");
            }
        }
    }
}

public int getIndexRuleWithBestFitness(HashMap<Integer, Double> rules) {
    double temp, bestValue = 0.0;
    int indexBestValue = -1;

    for (Integer key : rules.keySet()) {
        temp = rules.get(key);

        if (temp > bestValue) {
            bestValue = temp;
            indexBestValue = key;
        }
    }

    return indexBestValue;
}

public ArrayList<Double> getAUPRC(ArrayList<Integer> numberCoveredExamples, ArrayList<double[]> meanClassLabelVectorCovered, ArrayList<ArrayList<Integer>> indexCoveredExamples) {
    ArrayList<Double> AUPRC = new ArrayList<>();

    for (int i = 0; i < numberCoveredExamples.size(); i++) {
        //Obtain the predictions
        double[][] matrixPredictions = new double[numberCoveredExamples.get(i)][classes.length];

        for (int j = 0; j < numberCoveredExamples.get(i); j++) {
            System.arraycopy(meanClassLabelVectorCovered.get(i), 0, matrixPredictions[j], 0, classes.length);
        }

        //min Mips = 5
        if (indexCoveredExamples.get(i).size() >= 5) {
            AUPRC.add(i, evaluationAUPRCFitness(matrixPredictions, indexCoveredExamples.get(i)));
        }
    }

    return AUPRC;
}

public double evaluationAUPRCFitness(double[][] matrixPredictions, ArrayList<Integer> indexExamples) {

    //Store precision and recall values
    ArrayList<double[]> valuesPrecisionRecall = new ArrayList<double[]>();
    ArrayList<Double> thresholdValues = new ArrayList<Double>(Arrays.asList(0.0, 2.0, 4.0,
            6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0,
            46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0, 64.0, 66.0, 68.0, 70.0, 72.0, 74.0, 76.0, 78.0, 80.0, 82.0, 84.0,
            86.0, 88.0, 90.0, 92.0, 94.0, 96.0, 98.0, 100.0));

    double AUPRCFitness = 0;

    //Iterate over all thresholds
    for (int indexThres = 0; indexThres < thresholdValues.size(); indexThres++) {

        //Matrix to store the outputs on the test data after applying thresholds
        int[][] binaryMatrix = new int[indexExamples.size()][classes.length];

        //Threshold values used
        double threshold = thresholdValues.get(indexThres) / 100;

        //Apply the threshold
//            applyThresholds(binaryMatrix, matrixPredictions, threshold, 1);
        Evaluation.applyThresholds(binaryMatrix, matrixPredictions, threshold, 0);

        ArrayList<int[]> trueClasses = new ArrayList<int[]>();
        int[] classesB = new int[classes.length];

        for (int i = 0; i < indexExamples.size(); i++) {
            int posExample = indexExamples.get(i);
            System.arraycopy(binaryClassesTrain.get(posExample), 0, classesB, 0, classes.length);
            trueClasses.add(classesB);
        }

        //Hierarchical Precision and Recall evaluation metrics
        double[] evalResults = evaluationPrecRec(trueClasses, binaryMatrix);
        valuesPrecisionRecall.add(evalResults);

    }

    //Calculate AU(PRC)
    AUPRCFitness = calculateAUPRCFitness(valuesPrecisionRecall);

    return AUPRCFitness;

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
        ArrayList<ArrayList<Double>> points = Evaluation.getPoints(dataInterpolation, count);

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

public static double calculateAreaUnderCurve(ArrayList<ArrayList<Double>> recall,
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

public double[] evaluationPrecRec(ArrayList<int[]> trueClasses, int[][] predictedClasses) {
    //Store the results
    double[] evalResults = new double[5];

    //Sum of predicted and real classes
    double sumIntersection = 0;
    double sumPredicted = 0;
    double sumReal = 0;
    double FP = 0;

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

public ArrayList<Double> getVarianceGain(ArrayList<ArrayList<Integer>> indexCoveredExamples, ArrayList<ArrayList<Integer>> indexUncoveredExamples,
        ArrayList<Integer> numberCoveredExamples, ArrayList<double[]> meanClassLabelVectorCovered, ArrayList<double[]> meanClassLabelVectorUncovered) {

    ArrayList<Double> varianceGain = new ArrayList<>();

    for (int i = 0; i < indexCoveredExamples.size(); i++) {
        double[][] matrixPredictions = new double[numberCoveredExamples.get(i)][classes.length];

        for (int j = 0; j < numberCoveredExamples.get(i); j++) {
            System.arraycopy(meanClassLabelVectorCovered.get(i), 0, matrixPredictions[j], 0, classes.length);
        }

        //Variance gain of the set of all training examples
        double[] meanClassesLabelAll = meanClassLabelVectorAllClassesTrain;
        /*if (Parameters.getMultiLabel() == 0) {
             meanClassesLabelAll = Results.getHigherProbabilities(meanClassesLabelAll);
             }*/
        double varianceGainAll = getVarianceGain(datasetTrain.size(), meanClassesLabelAll);
        //Variance gain of the set of covered examples
        double varianceGainCovered = getVarianceGain(indexCoveredExamples.get(i), meanClassLabelVectorCovered.get(i));
        //Variance gain of the set of uncovered examples
        double varianceGainUncovered = getVarianceGain(indexUncoveredExamples.get(i), meanClassLabelVectorUncovered.get(i));

        int numTotalExamples = datasetTrain.size();
        double term11 = (double) indexCoveredExamples.get(i).size() / (double) numTotalExamples;
        double term1 = term11 * varianceGainCovered;
        double term22 = (double) indexUncoveredExamples.get(i).size() / (double) numTotalExamples;
        double term2 = term22 * varianceGainUncovered;

//        Variance Gain Fitness
//        -------------------------------------------------------------------------------
        double varianceGainTotal = varianceGainAll - term1 - term2;
//        double percentageCoverage = (double) numberCoveredExamples.get(i) / datasetTrain.size();
//        (0.4 * varianceGainTotal) + (0.6 * percentageCoverage) / percentageCoverage;
        varianceGain.add(i, varianceGainTotal);
    }

    return varianceGain;
}

public double getVarianceGain(int numExamples, double[] meanClassLabel) {
    double varianceGain = 0;
    double[] weights = weightingScheme;
    ArrayList<int[]> binaryClasses = binaryClassesTrain;

    for (int i = 0; i < numExamples; i++) {
        double sum = 0;
        for (int j = 0; j < meanClassLabel.length; j++) {
            sum += weights[j] * Math.pow((binaryClasses.get(i)[j] - meanClassLabel[j]), 2);
        }
        varianceGain += sum;
    }

    varianceGain = varianceGain / numExamples;

    return varianceGain;
}

public double getVarianceGain(ArrayList<Integer> indexExamples, double[] meanClassLabel) {

    double varianceGain = 0;
    double[] weights = weightingScheme;
    ArrayList<int[]> binaryClasses = binaryClassesTrain;

    for (int i = 0; i < indexExamples.size(); i++) {
        double sum = 0;
        int[] binaryVector = binaryClasses.get(indexExamples.get(i));
        for (int j = 0; j < meanClassLabel.length; j++) {
            sum += weights[j] * Math.pow((binaryVector[j] - meanClassLabel[j]), 2);
        }
        varianceGain += sum;
    }

    if (indexExamples.size() > 0) {
        varianceGain = varianceGain / indexExamples.size();
    }

    return varianceGain;
}

public static void main(String[] args) {
//    ClassificationWithRule c = new ClassificationWithRule();
//    c.run("makePredictions");
//    c.run("makePredictionsImproved");
//auxMinCov / generationReboots * 2
int result = 20 / (1 * 2);
        System.out.println("result: "+ result);
}
}
