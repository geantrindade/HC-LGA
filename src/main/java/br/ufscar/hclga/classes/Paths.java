package br.ufscar.hclga.classes;

import java.io.File;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Paths {

private static String nameFilePredictions = "predictions.txt";
private static String nameFileRules = "rules.txt";
private static String nameFilePRclasses = "PRclasses.txt";
private static String nameFileInterpolationClasses = "interpolationClasses.txt";
private static String nameFileInterpolation = "interpolation.txt";
private static String nameFileAUPRCclasses = "AUPRCclasses.txt";
private static String nameFileFmeasureLevels = "fmeasureLevels.txt";
private static String nameFileMeanFmeasureLevels = "meanFmeasureLevels.txt";
private static String nameFileMeanAUPRCclasses = "meanAUPRCclasses.txt";
private static String nameFileMeanAUPRC = "meanAUPRC.txt";
private static String nameFileTrainingTimes = "trainingTimes.txt";
private static String nameFileMeanTrainingTimes = "meanTrainingTimes.txt";

//Path where the datasets are storaged
private static String pathDatasets;

//A Results folder will be created in this path
private static String pathResults = "Results/";

/**
 * @param aNameFileTrainingTimes the nameFileTrainingTimes to set
 */
public static void setNameFileTrainingTimes(String aNameFileTrainingTimes) {
    nameFileTrainingTimes = aNameFileTrainingTimes;
}

public Paths() {
    pathDatasets = Parameters.getPathDatasets();
    pathResults = pathResults + "SingleLabel/";

    //Create directories to save results
    String pathResultsRuns = pathResults + Parameters.getHierarchyType() + "/" + Parameters.getFileDatasetTest() + "/";

    for (int i = 1; i <= Parameters.getNumberRuns(); i++) {
        new File(pathResultsRuns + "Run" + i).mkdirs();
    }
}

public static String getPathDatasets() {
    return pathDatasets;
}

public static String getPathResults() {
    return pathResults;
}

public static String getNameFileRules() {
    return nameFileRules;
}

public static String getNameFilePredictions() {
    return nameFilePredictions;
}

public static String getNameFilePRclasses() {
    return nameFilePRclasses;
}

public static String getNameFileInterpolationClasses() {
    return nameFileInterpolationClasses;
}

public static String getNameFileInterpolation() {
    return nameFileInterpolation;
}

public static String getNameFileAUPRCclasses() {
    return nameFileAUPRCclasses;
}

public static String getNameFileMeanAUPRCclasses() {
    return nameFileMeanAUPRCclasses;
}

public static String getNameFileMeanAUPRC() {
    return nameFileMeanAUPRC;
}

public static String getNameFileMeanFmeasureLevels() {
    return nameFileMeanFmeasureLevels;
}

public static String getNameFileFmeasureLevels() {
    return nameFileFmeasureLevels;
}

public static String getNameFileTrainingTimes() {
    return nameFileTrainingTimes;
}

public static String getNameFileMeanTrainingTimes() {
    return nameFileMeanTrainingTimes;
}

}
