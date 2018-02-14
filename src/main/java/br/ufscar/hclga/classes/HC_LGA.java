package br.ufscar.hclga.classes;

import java.util.ArrayList;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class HC_LGA {

public static int numRun;
private static ArrayList<double[]> meansAllAUPRCClasses;
private static double[] meanAUPRCTest;
private static ArrayList<double[]> meanFmeasures;

/**
 * @param args the command line arguments
 */
public static void main(String[] args) {
    //Get the parameters
//    new Parameters("C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\configFile.txt");
    new Parameters("configFile.txt");

    //Create directories for results
    new Paths();

    GeneticOperators.multiObj = false;

    //Store the AUPRCs for each run
    ArrayList<ArrayList<double[]>> AUPRCsRunsClasses = new ArrayList<ArrayList<double[]>>();
    double[] allAUPRCTest = new double[Parameters.getNumberRuns()];

    //Store the fmeasures for each run
    ArrayList<double[]> fmeasureRuns = new ArrayList<double[]>();

    //Will save the probability of using clausule because it can be modified during execution
    double probabilityUseClausule = Parameters.getProbabilityUseClausule();

    //Will store size of the test dataset to calculate means
    int sizeDatasetTest = 0;

    //Store the training times of all executions
    double[] allMsTimes = new double[Parameters.getNumberRuns()];
    double[] allSTimes = new double[Parameters.getNumberRuns()];
    double[] allMTimes = new double[Parameters.getNumberRuns()];
    double[] allHTimes = new double[Parameters.getNumberRuns()];

    //Execute the algorithm as many times as the user wants
    for (numRun = 1; numRun <= Parameters.getNumberRuns(); numRun++) {
        //Will save training times
        double[] trainingTimes = new double[4];

        //Evaluate the predictions in the test dataset
        Evaluation evaluation = new Evaluation();

        // ================================ begin Obj1 ===================================
        //Set the probability of using clausule in the beggining of each executuon
        Parameters.setProbabilityUseClausule(probabilityUseClausule);

        //Read the training and validation datasets
        Datasets datasets = new Datasets();

        //Mark training time
        Chronometer chron = new Chronometer();
        chron.start();

        //Build the binary structure to store the classes
        Classes.buildClassesStructureTrain();

        double[] defaultRule = Classes.getMeanClassLabelVectorAllClasses().clone();

        if (Parameters.getMultiLabel() == 0) {
            defaultRule = Results.getHigherProbabilities(defaultRule);
        }

//        int numberAttempt = 0;
//        int maxAttempts = 0;

        ArrayList<Individual> bestIndividualsObj1 = new ArrayList<Individual>();

        int maxUncoveredExamples = Parameters.getMaxUncoveredExamples();

        while (Datasets.getDatasetTrain().size() > maxUncoveredExamples) {
//            numberAttempt++;
//            maxAttempts++;

            Population populationObj1 = new Population(false);

//            System.out.println("rules obj1 created");
            //Initiate evolution          
            Evolution evolution = new Evolution(populationObj1, false);

//            System.out.println("evolution " + numberAttempt);
            //Gets best rules
            Individual bestIndividual = evolution.getBestIndividual();

            //Verify if the best rule covers a minimum specified number of examples
//            if (bestIndividual.getNumberCoveredExamples() >= Parameters.getMinCoveredExamplesRule()) {
                //&& bestIndividual.getNumberCoveredExamples() <= Parameters.getMaxCoveredExamplesRule()) {

                //if (Parameters.getMultiLabel() == 1) {
            bestIndividualsObj1.add(bestIndividual);

                //Remove covered examples from the dataset
            ArrayList<Integer> indexesCoveredExamples = bestIndividual.getIndexCoveredExamples();

            int trainSize = Datasets.getDatasetTrain().size();

            System.out.println("hF: " + bestIndividual.getHFmeasure());
            
            Datasets.removeTrainExamples(indexesCoveredExamples);
             
            System.out.println("number of train examples covered: " + (trainSize - Datasets.getDatasetTrain().size()));
             
            //Build the binary structure to store the classes of the training data
             Classes.buildClassesStructureTrain();

            //maxAttempts = 0;
            
            
            //System.out.println("current number of initial rules: " + Parameters.getNumberInitialRules());
            Double newValue = 0.0; 
            
            if(Parameters.getDatasetTrain().contains("mips")){
                newValue = Datasets.getDatasetTrain().size() * 0.03;
            }else if(Parameters.getDatasetTrain().contains("repbase")){
                newValue = Datasets.getDatasetTrain().size() * 0.02;
            }
             
            if (newValue < 100) { //min number of rules
                Parameters.setNumberInitialRules(100);
            } else {
                Parameters.setNumberInitialRules(newValue.intValue());
            }
//            }

//                System.out.println("new number of initial rules: " + Parameters.getNumberInitialRules());
            //If it stays X attempts without removing examples, let's
            //decrease the number of initial terms by 2
//            int approxNumUsedTerms = (int) (Datasets.getInfoAttributes().size() * Parameters.getProbabilityUseClausule());
//            int newApproxNumberUsedTerms = approxNumUsedTerms - 2;
////                if (maxAttempts == 50 && newApproxNumberUsedTerms > 0) {
//            if (maxAttempts == Parameters.getNumberGenerations() / 2 && newApproxNumberUsedTerms > 0) {
//                double newPercentage = (double) (newApproxNumberUsedTerms * 100) / (Datasets.getInfoAttributes().size() * 100);
//                Parameters.setProbabilityUseClausule(newPercentage);
//                maxAttempts = 0;
//            }

            //let the remaining examples uncovered
//            if (maxAttempts == Parameters.getNumberGenerations()) {
//                maxUncoveredExamples = Datasets.getDatasetTrain().size();
//                System.out.println("para de procurar regras");
//            }

            populationObj1 = null;
            evolution = null;
            bestIndividual = null;
            System.gc();
        }
// ================================ end Obj1 ===================================

        System.out.println("\n============= fim obj1 ============");
        chron.stop();

        System.out.println("segundos: " + chron.stime());
        System.out.println("minutos: " + chron.mtime());
        System.out.println("horas: " + chron.htime());
        System.out.println("\nnumero de regras: " + bestIndividualsObj1.size());
        System.out.println("============= fim obj1 ============");

        trainingTimes[0] = chron.time();
        trainingTimes[1] = chron.stime();
        trainingTimes[2] = chron.mtime();
        trainingTimes[3] = chron.htime();

        String bk = Paths.getNameFileTrainingTimes();

        Paths.setNameFileTrainingTimes("objective1.txt");
        Results.saveTrainingTimesRun(trainingTimes);

        Paths.setNameFileTrainingTimes(bk);

        chron.resume();

//        Results.saveRules(bestIndividualsObj1, "C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\obj1\\");
        Results.saveRules(bestIndividualsObj1, "obj1/");

//        ClassificationWithRule cl = new ClassificationWithRule("mips", 1, 10,
//                "C:\\Users\\gean_\\Desktop\\", "C:\\Users\\gean_\\Desktop\\", "C:\\Users\\gean_\\Documents"
//                + "\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\obj1\\",
//                "C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\obj1\\");
//
//        cl.run("makePredictions");




        //================================================================================
        // ================================ begin Obj2 ===================================
        //================================================================================
        //Disable mutation because this can make the search go to other places outside of the interesting region already found
        Parameters.setMutationRate(0.0);
        Parameters.setProbabilityUseClausule(0.0);
        
        Parameters.setMinFitnessTh(0.8); //higher minTh
        Parameters.setMaxFitnessTh(0.9); //higher maxTh
//        Parameters.setNumberAttempts(Parameters.getNumberAttempts() * 2);
//        Parameters.setNumberGenerations(Parameters.getNumberGenerations() * 2);
        
        GeneticOperators.multiObj = true;

        //Read the training AGAIN
        datasets = new Datasets();

        //Build the binary structure to store the classes AGAIN
        Classes.buildClassesStructureTrain();

        defaultRule = Classes.getMeanClassLabelVectorAllClasses().clone();

        if (Parameters.getMultiLabel() == 0) {
            defaultRule = Results.getHigherProbabilities(defaultRule);
        }

//        maxAttempts = 0;

        ArrayList<Individual> bestIndividualsObj2 = new ArrayList<Individual>();

        maxUncoveredExamples = Parameters.getMaxUncoveredExamples();

        while (Datasets.getDatasetTrain().size() > maxUncoveredExamples) {
//            maxAttempts++;

            Population populationObj2 = new Population(true);

            ArrayList<double[]> popObj1 = new ArrayList<>();
            ArrayList<ArrayList<Integer>> activeTermsPopObj1 = new ArrayList<>();

            for (int i = 0; i < bestIndividualsObj1.size(); i++) {
                popObj1.add(bestIndividualsObj1.get(i).getRule());
                activeTermsPopObj1.add(bestIndividualsObj1.get(i).getPosActiveTerms());
            }

            populationObj2.setPopulation(popObj1);
            populationObj2.setActiveTerms(activeTermsPopObj1);

//            System.out.println("rules obj2 created");
            //Initiate evolution          
            Evolution evolution = new Evolution(populationObj2, true);

//            System.out.println("evolution " + numberAttempt);
            //Gets best rules
            Individual bestIndividual = evolution.getBestIndividual();

            //Verify if the best rule covers a minimum specified number of examples
//            if (bestIndividual.getNumberCoveredExamples() >= Parameters.getMinCoveredExamplesRule()) {
            bestIndividualsObj2.add(bestIndividual);

                //Remove covered examples from the dataset
            ArrayList<Integer> indexesCoveredExamples = bestIndividual.getIndexCoveredExamples();

            int trainSize = Datasets.getDatasetTrain().size();

            System.out.println("hF: " + bestIndividual.getHFmeasure());
            
            Datasets.removeTrainExamples(indexesCoveredExamples);
            
            System.out.println("number of train examples covered: " + (trainSize - Datasets.getDatasetTrain().size()));
            
            Classes.buildClassesStructureTrain();

//                maxAttempts = 0;
//            }

            //If it stays X attempts without removing examples, let's
            //decrease the number of initial terms by 2
//            int approxNumUsedTerms = (int) (Datasets.getInfoAttributes().size() * Parameters.getProbabilityUseClausule());
//            int newApproxNumberUsedTerms = approxNumUsedTerms - 2;
//                if (maxAttempts == 50 && newApproxNumberUsedTerms > 0) {
//            if (maxAttempts == Parameters.getNumberGenerations() / 2 && newApproxNumberUsedTerms > 0) {
//                double newPercentage = (double) (newApproxNumberUsedTerms * 100) / (Datasets.getInfoAttributes().size() * 100);
//                Parameters.setProbabilityUseClausule(newPercentage);
//                maxAttempts = 0;
//            }

            //After 100 executions without finding a rule that covers the examples,
            //let the remaining examples uncovered
//            if (maxAttempts == Parameters.getNumberGenerations()) {
//                maxUncoveredExamples = Datasets.getDatasetTrain().size();
////                System.out.println("para de procurar regras");
//            }

            populationObj2 = null;
            evolution = null;
            bestIndividual = null;
            System.gc();
        }

// ================================ end Obj2 ===================================
//Save training time
        chron.stop();

        System.out.println("\n============= fim obj2 ============");
        System.out.println("segundos: " + chron.stime());
        System.out.println("minutos: " + chron.mtime());
        System.out.println("horas: " + chron.htime());
        System.out.println("\nnumero de regras: " + bestIndividualsObj2.size());
        System.out.println("============= fim obj2 ============");

//        Results.saveRules(bestIndividualsObj2, "C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\obj2\\");
        Results.saveRules(bestIndividualsObj2, "obj2/");

//        ClassificationWithRule cl2 = new ClassificationWithRule("mips", 1, 10,
//                "C:\\Users\\gean_\\Desktop\\", "C:\\Users\\gean_\\Desktop\\", "C:\\Users\\gean_\\Documents"
//                + "\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\obj2\\",
//                "C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\obj2\\");
        trainingTimes[0] = chron.time();
        trainingTimes[1] = chron.stime();
        trainingTimes[2] = chron.mtime();
        trainingTimes[3] = chron.htime();

        bk = Paths.getNameFileTrainingTimes();

        Paths.setNameFileTrainingTimes("objective2.txt");
        Results.saveTrainingTimesRun(trainingTimes);

        Paths.setNameFileTrainingTimes(bk);

        Results.saveTrainingTimesRun(trainingTimes);

//        Results.saveTrainingTimesRun(trainingTimes);
//        allMsTimes[numRun - 1] = trainingTimes[0];
//        allSTimes[numRun - 1] = trainingTimes[1];
//        allMTimes[numRun - 1] = trainingTimes[2];
//        allHTimes[numRun - 1] = trainingTimes[3];
//Load test data
//        datasets.readTestData(Paths.getPathDatasets() + Parameters.getFileDatasetTest());
//        Classes.buildClassesStructureTest();
//        sizeDatasetTest = Datasets.getDatasetTest().size();
/*if(Datasets.getDatasetTest().size() == 0 || bestIndividuals.size() == 0){
                System.out.println();
            }*/
//Order the indivuduals by the fitness in the whole training set
//---------------------------------------------------------------------------------------------------------------------
/*new Datasets();
             Classes.buildClassesStructureTrain();

             ArrayList<Individual> finalIndividuals = new ArrayList<Individual>();
             for (int i = 0; i < bestIndividuals.size(); i++) {
             finalIndividuals.add(new Individual(bestIndividuals.get(i).getRule(), bestIndividuals.get(i).getPosActiveTerms()));
             }

             Collections.sort(finalIndividuals);*/
//---------------------------------------------------------------------------------------------------------------------
//Obtain the predictions
//        double[][] matrixPredictions = Results.obtainPredictions(bestIndividualsObj2, defaultRule);
//        //For single-label, lets get only the higher probabilities in the vectors for each level
//        if (Parameters.getMultiLabel() == 0) {
//            matrixPredictions = Results.getHigherProbabilities(matrixPredictions);
//
//            //Lets calculate other evaluation measures for single-label problems
//            //F-measure per level
//            evaluation.evaluationFmeasure(matrixPredictions);
//            double[] fmeasureLevels = evaluation.getFmeasureLevels();
//
//            //Save the fmeasures for this run
//            Results.saveFmeasureRun(fmeasureLevels);
//            fmeasureRuns.add(fmeasureLevels);
//        }
//
//        evaluation.evaluationAUPRC(matrixPredictions);
//        double AUPRC = evaluation.getAUPRC();
//        ArrayList<double[]> AUPRCClasses = evaluation.getAUPRCClasses();
//
//        //Save predictions, rules and AUPRC values
//        Results.savePredictions(matrixPredictions, AUPRC, AUPRCClasses);
//
//        String predictions = "";
//        for (int i = 0; i < matrixPredictions.length; i++) {
//            for (int j = 0; j < matrixPredictions[i].length; j++) {
//                predictions += matrixPredictions[i][j] + " ";
//            }
//            predictions += "\n";
//        }
//
//        String realTestClasses = "";
//        int[] aux;
//        for (int i = 0; i < Classes.getBinaryClassesTest().size(); i++) {
//            aux = Classes.getBinaryClassesTest().get(i);
//
//            for (int j = 0; j < aux.length; j++) {
//                realTestClasses += aux[j] + " ";
//            }
//            realTestClasses += "\n";
//        }
//        Validacao v = new Validacao();
//        String eval = v.run(realTestClasses, predictions, 0.5);
//        System.out.println("eval: " + eval);
//        v.escreveArquivo("C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\evaluation.txt", eval);
//        v.escreveArquivo("evaluation.txt", eval);
//        Results.saveRules(bestIndividualsObj2);
//        allAUPRCTest[numRun - 1] = AUPRC;
//        AUPRCsRunsClasses.add(AUPRCClasses);
//Free test dataset
//        datasets.freeTestDataset();
    }

//Calculate the meanAUPRC for each class and overall
//    meansAllAUPRCClasses = Results.calculateMeansAUPRCClasses(AUPRCsRunsClasses);
//    meanAUPRCTest = Results.calculateMeanSd(allAUPRCTest);
//Calculate mean execution times
//    double[] meanMsTimes = Results.calculateMeanSd(allMsTimes);
//    double[] meanSTimes = Results.calculateMeanSd(allSTimes);
//    double[] meanMTimes = Results.calculateMeanSd(allMTimes);
//    double[] meanHTimes = Results.calculateMeanSd(allHTimes);
//Save mean AUPRC for each class and overall
//    Results.saveMeansAUPRCClasses(meansAllAUPRCClasses, meanAUPRCTest, allAUPRCTest, sizeDatasetTest);
//Save mean training times
//    Results.saveMeanTrainingTimes(meanMsTimes, meanSTimes, meanMTimes, meanHTimes);
//    if (Parameters.getMultiLabel()
//            == 0) {
//        //Calculate mean and sd for all executions in each level
//        meanFmeasures = Results.calculateMeanSdFmesureLevels(fmeasureRuns);
//        //Save the mean Fmeasure for each level
//        Results.saveMeanFmeasureLevels(meanFmeasures);
//    }
//        Validacao v = new Validacao();
//        v.run();
}

public static int getNumRun() {
    return numRun;
}
}
