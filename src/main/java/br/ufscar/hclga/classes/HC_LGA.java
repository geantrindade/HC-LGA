package br.ufscar.hclga.classes;

import java.util.ArrayList;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class HC_LGA {

public static void main(String[] args) {
    //Get the parameters
//    new Parameters("C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\configFile.txt");
    new Parameters("configFile.txt");

    //Create directories for results
    new Paths();

    GeneticOperators.multiObj = false;

    //Will save the probability of using clausule because it can be modified during execution
    double probabilityUseClausule = Parameters.getProbabilityUseClausule();

    //Will store size of the test dataset to calculate means
    int sizeDatasetTest = 0;

    //Execute the algorithm as many times as the user wants
    for (int numRun = 1; numRun <= Parameters.getNumberRuns(); numRun++) {
        //Will save training times
        double[] trainingTimes = new double[4];

        //Evaluate the predictions in the test dataset
        Evaluation evaluation = new Evaluation();

        //================================================================================
        // ================================ begin Obj1 ===================================
        //================================================================================
        //
        //Set the probability of using clausule in the beggining of each executuon
        Parameters.setProbabilityUseClausule(probabilityUseClausule);

        //Read the training and validation datasets
        Datasets datasets = new Datasets();

        //Mark training time
        Chronometer chron = new Chronometer();
        chron.start();

        //Build the binary structure to store the classes
        Classes.buildClassesStructureTrain();

        //will store solutions for the first objective
        ArrayList<Individual> bestIndividualsObj1 = new ArrayList<>();

        int maxUncoveredExamples = Parameters.getMaxUncoveredExamples();

        while (Datasets.getDatasetTrain().size() > maxUncoveredExamples) {
            Population populationObj1 = new Population(false);

            Evolution evolution = new Evolution(populationObj1, false);

            //Get best rules
            Individual bestIndividual = evolution.getBestIndividual();
            bestIndividualsObj1.add(bestIndividual);

            //Remove covered examples from the dataset
            ArrayList<Integer> indexesCoveredExamples = bestIndividual.getIndexCoveredExamples();

            System.out.println("hF: " + bestIndividual.getHFmeasure());

            int trainSize = Datasets.getDatasetTrain().size();

            Datasets.removeTrainExamples(indexesCoveredExamples);

            System.out.println("number of train examples covered: " + (trainSize - Datasets.getDatasetTrain().size()));

            //Build the binary structure to store the classes of the training data
            Classes.buildClassesStructureTrain();

            Double newValue = 0.0;

            if (Parameters.getDatasetTrain().contains("mips")) {
                newValue = Datasets.getDatasetTrain().size() * 0.02;
            } else if (Parameters.getDatasetTrain().contains("repbase")) {
                newValue = Datasets.getDatasetTrain().size() * 0.02;
            }

            if (newValue < 100) { //min number of rules
                Parameters.setNumberInitialRules(100);
            } else {
                Parameters.setNumberInitialRules(newValue.intValue());
            }

            //clean cache
            populationObj1 = null;
            evolution = null;
            bestIndividual = null;
            System.gc();
        }
// ================================ end Obj1 ===================================

        System.out.println("\n============= end obj1 ============");
        chron.stop();

        System.out.println("seconds: " + chron.stime());
        System.out.println("minutes: " + chron.mtime());
        System.out.println("hours: " + chron.htime());
        System.out.println("\nnumber of rules: " + bestIndividualsObj1.size());
        System.out.println("============= end obj1 ============");

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
//
//
        //================================================================================
        // ================================ begin Obj2 ===================================
        //================================================================================
        //
        //
        //Disable mutation because this can make the search go to other places outside of the interesting region already found
        Parameters.setMutationRate(0.0);
        Parameters.setProbabilityUseClausule(0.0);

        Parameters.setMinFitnessTh(0.8); //higher minTh
        Parameters.setMaxFitnessTh(0.9); //higher maxTh
        Parameters.setNumberAttempts(Parameters.getNumberAttempts() * 2);
        Parameters.setNumberGenerations(Parameters.getNumberGenerations() * 2);

        GeneticOperators.multiObj = true;

        //Read the training AGAIN
        datasets = new Datasets();

        //Build the binary structure to store the classes AGAIN
        Classes.buildClassesStructureTrain();

        //will store solutions for the second objective
        ArrayList<Individual> bestIndividualsObj2 = new ArrayList<Individual>();

        maxUncoveredExamples = Parameters.getMaxUncoveredExamples();

        while (Datasets.getDatasetTrain().size() > maxUncoveredExamples) {
            Population populationObj2 = new Population(true);

            ArrayList<double[]> popObj1 = new ArrayList<>();
            ArrayList<ArrayList<Integer>> activeTermsPopObj1 = new ArrayList<>();

            for (int i = 0; i < bestIndividualsObj1.size(); i++) {
                popObj1.add(bestIndividualsObj1.get(i).getRule());
                activeTermsPopObj1.add(bestIndividualsObj1.get(i).getPosActiveTerms());
            }

            populationObj2.setPopulation(popObj1);
            populationObj2.setActiveTerms(activeTermsPopObj1);

            Evolution evolution = new Evolution(populationObj2, true);

            Individual bestIndividual = evolution.getBestIndividual();
            bestIndividualsObj2.add(bestIndividual);

            ArrayList<Integer> indexesCoveredExamples = bestIndividual.getIndexCoveredExamples();

            int trainSize = Datasets.getDatasetTrain().size();

            System.out.println("hF: " + bestIndividual.getHFmeasure());

            Datasets.removeTrainExamples(indexesCoveredExamples);

            System.out.println("number of train examples covered: " + (trainSize - Datasets.getDatasetTrain().size()));

            Classes.buildClassesStructureTrain();

            populationObj2 = null;
            evolution = null;
            bestIndividual = null;
            System.gc();
        }
// ================================ end Obj2 ===================================
//
        //Save training time
        chron.stop();

        System.out.println("\n============= end obj2 ============");
        System.out.println("seconds: " + chron.stime());
        System.out.println("minutes: " + chron.mtime());
        System.out.println("hours: " + chron.htime());
        System.out.println("\nnumber of rules: " + bestIndividualsObj2.size());
        System.out.println("============= end obj2 ============");

//        Results.saveRules(bestIndividualsObj2, "C:\\Users\\gean_\\Documents\\NetBeansProjects\\HC-LGA\\src\\main\\java\\br\\ufscar\\hclga\\obj2\\");
        Results.saveRules(bestIndividualsObj2, "obj2/");

        trainingTimes[0] = chron.time();
        trainingTimes[1] = chron.stime();
        trainingTimes[2] = chron.mtime();
        trainingTimes[3] = chron.htime();

        bk = Paths.getNameFileTrainingTimes();

        Paths.setNameFileTrainingTimes("objective2.txt");
        Results.saveTrainingTimesRun(trainingTimes);

        Paths.setNameFileTrainingTimes(bk);

        Results.saveTrainingTimesRun(trainingTimes);
    }
}
}
