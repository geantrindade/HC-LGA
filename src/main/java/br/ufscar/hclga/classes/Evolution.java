package br.ufscar.hclga.classes;

import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Evolution {

private ArrayList<Individual> currentPopulation;
private Individual bestIndividual;

public Evolution(Population population, boolean multiObj) {
    currentPopulation = new ArrayList<Individual>();
    ArrayList<Individual> backupInitialPopulation = new ArrayList<Individual>();

    //Initializing current population
    int populationSize = Parameters.getNumberInitialRules();
    int obj = 0;

    if (multiObj == false) {
        obj = 1;

        for (int i = 0; i < populationSize; i++) {
            currentPopulation.add(new Individual(population.getPopulation().get(i), population.getActiveTerms().get(i), multiObj));
        }

    } else {
        obj = 2;

        for (int i = 0; i < population.getPopulation().size(); i++) {
            currentPopulation.add(new Individual(population.getPopulation().get(i), population.getActiveTerms().get(i), multiObj));
        }
    }

    Collections.sort(currentPopulation);
    bestIndividual = new Individual(currentPopulation.get(0).getRule().clone(), currentPopulation.get(0).getPosActiveTerms(), multiObj);

    backupInitialPopulation = currentPopulation; //storage for future check
    ArrayList<Individual> nextPopulation = new ArrayList<Individual>();

    ArrayList<Individual> children;

    int numGenerations = 1;
    int attempts = 1;
    int generationReboots = 0;
    int paretoReboots = 0;

    //prevent useless iterations
    if (Parameters.getMinCoveredExamplesRule() > Datasets.getDatasetTrain().size()) {
        Parameters.setMinCoveredExamplesRule(Parameters.getMinCoveredExamplesRule() / 2);
    }

    // attempts = number of attempts to find a better rule that the current better rule
    System.out.println("\n\n=============== obj" + obj + " ===================");
    System.out.println("remaining examples: " + Datasets.getDatasetTrain().size());
    System.out.println("initial rules: " + currentPopulation.size());

    //storage the original values
    int auxMinCov = Parameters.getMinCoveredExamplesRule();
    int auxMaxCov = Parameters.getMaxCoveredExamplesRule();
    int auxTournament = Parameters.getSizeTournament();

    while (numGenerations <= Parameters.getNumberGenerations() && attempts <= Parameters.getNumberAttempts()) {
        System.out.println("\ngeneration " + numGenerations);
        System.out.println("attempt " + attempts);
        System.out.println("fitness best rule: " + bestIndividual.getFitness());

        // max value reached, so dont need to keeping searching
        if (bestIndividual.getHFmeasure() >= Parameters.getMaxFitnessTh()) {
            System.out.println("\nmax value reached!!!");

            for (Individual i : currentPopulation) { //check if other individuals got good values
                if (i.getHFmeasure() >= Parameters.getMaxFitnessTh()) {
                    //for objective 1
                    if (multiObj == false) {

                        if (i.getFitness(0) >= bestIndividual.getFitness(0)) {

                            if (i.getFitness(1) > bestIndividual.getFitness(1)) {
                                System.out.println("similar or better value for obj1!");
                                System.out.println("HIGHER value for obj2!");

                                System.out.println("old bestRule obj1: " + bestIndividual.getFitness(0));
                                System.out.println("old bestRule obj2: " + bestIndividual.getFitness(1));
                                System.out.println("new bestRule obj1: " + i.getFitness(0));
                                System.out.println("new bestRule obj2: " + i.getFitness(1));

                                bestIndividual = i;

                            } else if ((i.getFitness(1) == bestIndividual.getFitness(1))
                                    && (i.getNumberCoveredExamples() > bestIndividual.getNumberCoveredExamples())) {
                                System.out.println("similar or better value for obj1!");
                                System.out.println("similar value for obj2, but covered MORE examples!");

                                System.out.println("old bestRule n covered: " + bestIndividual.getNumberCoveredExamples());
                                System.out.println("new bestRule n covered: " + i.getNumberCoveredExamples());

                                bestIndividual = i;
                            }
                        }

                        //for objective 2
                    } else if (i.getFitness(1) >= bestIndividual.getFitness(1)) {
                        if (i.getFitness(0) > bestIndividual.getFitness(0)) {
                            System.out.println("similar or better value for obj2!");
                            System.out.println("HIGHER value for obj1!");

                            System.out.println("old bestRule obj1: " + bestIndividual.getFitness(0));
                            System.out.println("old bestRule obj2: " + bestIndividual.getFitness(1));
                            System.out.println("new bestRule obj1: " + i.getFitness(0));
                            System.out.println("new bestRule obj2: " + i.getFitness(1));

                            bestIndividual = i;

                        } else if ((i.getFitness(0) == bestIndividual.getFitness(0))
                                && (i.getNumberCoveredExamples() > bestIndividual.getNumberCoveredExamples())) {
                            System.out.println("similar or better value for obj2!");
                            System.out.println("similar value for obj1, but covered MORE examples!");

                            System.out.println("old bestRule n covered: " + bestIndividual.getNumberCoveredExamples());
                            System.out.println("new bestRule n covered: " + i.getNumberCoveredExamples());

                            bestIndividual = i;
                        }
                    }

                } else {
                    break;
                }
            }

            //dont have to iterate all the current population
            numGenerations = Parameters.getNumberGenerations();

            //keep going with the search
        } else {
            System.out.println("\nevolution started...");

            nextPopulation.addAll(GeneticOperators.elitism(this));
            children = GeneticOperators.uniformCrossoverDistance(this);
            children = GeneticOperators.mutation(children);
            nextPopulation.addAll(children);

            nextPopulation = GeneticOperators.localSearchOperatorMinMax(nextPopulation);

            currentPopulation.clear();
            currentPopulation.addAll(nextPopulation);
            nextPopulation.clear();
            Collections.sort(currentPopulation);

            System.out.println("\noldBestRule: " + bestIndividual.getFitness());
            System.out.println("newBestRule: " + currentPopulation.get(0).getFitness());

            if (currentPopulation.get(0).getFitness() > bestIndividual.getFitness()) {
                bestIndividual = null;
                System.gc();
                bestIndividual = currentPopulation.get(0);

                attempts = 0; //reset attempts
                System.out.println("\nnew bestRule found!!!");
                System.out.println("fitness: " + bestIndividual.getFitness());
            }
        }

        //last generation
        if (numGenerations == Parameters.getNumberGenerations()) {
            if (bestIndividual.getHFmeasure() < Parameters.getMinFitnessTh() && generationReboots < 3) {
                System.out.println("\nmin fitness threshold not achieved!!!");

                currentPopulation = backupInitialPopulation;
                Collections.sort(currentPopulation);

                //put the current bestRule into the new population in place to some worst rule
                currentPopulation.set(currentPopulation.size() - 5, bestIndividual);

                //update the current bestRule
                bestIndividual = new Individual(currentPopulation.get(0).getRule().clone(), currentPopulation.get(0).getPosActiveTerms(), multiObj);
                generationReboots++;

                System.out.println("restarting the evolutionary process with the initial population but "
                        + "with new parameters...");

                //reestart the evolution with new initial values
                if (Parameters.getMinCoveredExamplesRule() > 1
                        && (auxMinCov / (generationReboots * 2)) >= 1) {
                    Parameters.setMinCoveredExamplesRule(auxMinCov / (generationReboots * 2));
                    System.out.println("minCov decreased!");

                } else if (Parameters.getMaxCoveredExamplesRule() > 10
                        && (auxMaxCov / (generationReboots * 2)) >= 10) {
                    System.out.println("minCov cant be decreased anymore!!!");
                    Parameters.setMaxCoveredExamplesRule(auxMaxCov / (generationReboots * 2));
                    System.out.println("maxCov decreased!");

                    if (Parameters.getSizeTournament() > 2) {
                        Parameters.setSizeTournament(Parameters.getSizeTournament() - 1);
                        System.out.println("tournament decreased!");
                    }

                } else {
                    System.out.println("minCov, maxCov AND tournament size cant be decreased anymore!!!");
                    Parameters.setMinCoveredExamplesRule(1);
                    Parameters.setMaxCoveredExamplesRule(10);
                    Parameters.setSizeTournament(2);
                    System.out.println("setting them to the minimum...");

                }

                //reset 
                numGenerations = 0;
                attempts = 0;

                System.out.println("generation reboot " + generationReboots);

            } else if (bestIndividual.getHFmeasure() < Parameters.getMinFitnessTh() && generationReboots == 3) {
                System.out.println("\nnumber of generation reboots exhausted... evolution will be finished");

                //skip next checks
                numGenerations = Parameters.getNumberGenerations();
                attempts = Parameters.getNumberAttempts() + 1;

            } else {
                generationReboots = 0;
                System.out.println("\nmin fitness threshold achieved!!!");
            }
        }

        //last attempt
        if (attempts == Parameters.getNumberAttempts()) {
            if (bestIndividual.getHFmeasure() < Parameters.getMinFitnessTh()) {
                System.out.println("\nmin fitness threshold not achieved!!!");

                //decrease the minCovExamplesPerRule parameter        
                if (Parameters.getMinCoveredExamplesRule() > 5) {
                    Parameters.setMinCoveredExamplesRule(Parameters.getMinCoveredExamplesRule() - 5);
                    System.out.println("minCov decreased!");

                } else if (Parameters.getMaxCoveredExamplesRule() > 10) {
                    System.out.println("minCov cant be decreased anymore!!!");
                    Parameters.setMaxCoveredExamplesRule(Parameters.getMaxCoveredExamplesRule() / 2);
                    System.out.println("maxCov decreased!");

                    if (Parameters.getSizeTournament() > 2) {
                        Parameters.setSizeTournament(Parameters.getSizeTournament() - 1);
                        System.out.println("tournament decreased!");
                    }

                } else {
                    System.out.println("minCov, maxCov AND tournament size cant be decreased anymore!!!");
                    Parameters.setMinCoveredExamplesRule(1);
                    Parameters.setMaxCoveredExamplesRule(10);
                    Parameters.setSizeTournament(2);
                    System.out.println("setting them to the minimum...");
                }

                //reset
                attempts = 0;

            } else {
                System.out.println("\nmin fitness threshold achieved!!!");
            }
        }

        if ((numGenerations == Parameters.getNumberGenerations()
                || attempts == Parameters.getNumberAttempts())
                && multiObj == true) {
            //check if all the solutions in obj2 are in the pareto front
            boolean pareto = true;

            for (int i = 0; i < backupInitialPopulation.size(); i++) {
                if (bestIndividual.getFitness(0) < backupInitialPopulation.get(i).getFitness(0)
                        || bestIndividual.getFitness(1) < backupInitialPopulation.get(i).getFitness(1)) {

                    if (paretoReboots < 100) {
                        System.out.println("\nsolution not in the PF, rebooting evolution...");
                        paretoReboots++;
                        System.out.println("pareto reboot " + paretoReboots);
                        pareto = false;

                        currentPopulation = backupInitialPopulation; //reset  
                        Collections.sort(currentPopulation);

                        //put the current bestRule into the new population in place to some worst rule
                        currentPopulation.set(currentPopulation.size() - 5, bestIndividual);

                        bestIndividual.setFitness(0.0, 1); //force the evolution

                        if (Parameters.getMinCoveredExamplesRule() > 1
                                && (auxMinCov / (paretoReboots * 2)) >= 1) {
                            Parameters.setMinCoveredExamplesRule(auxMinCov / (paretoReboots * 2));
                            System.out.println("minCov decreased!");

                        } else if (Parameters.getMaxCoveredExamplesRule() > 10
                                && (auxMaxCov / (paretoReboots * 2)) >= 10) {
                            System.out.println("minCov cant be decreased anymore!!!");
                            Parameters.setMaxCoveredExamplesRule(auxMaxCov / (paretoReboots * 2));
                            System.out.println("maxCov decreased!");

                            if (Parameters.getSizeTournament() > 2) {
                                Parameters.setSizeTournament(Parameters.getSizeTournament() - 1);
                                System.out.println("tournament decreased!");
                            }

                        } else {
                            System.out.println("minCov, maxCov AND tournament size cant be decreased anymore!!!");
                            Parameters.setMinCoveredExamplesRule(1);
                            Parameters.setMaxCoveredExamplesRule(10);
                            Parameters.setSizeTournament(2);
                            System.out.println("setting them to the minimum...");
                        }

                        numGenerations = 0;
                        attempts = 0;

                        break;

                    } else {
                        System.out.println("\nnumber of reboots was exceeded, picking the best rule found...");
                        numGenerations = Parameters.getNumberGenerations();

                        break;
                    }
                }
            }

            if (pareto) {
                System.out.println("\nsolution in the pareto front!!!");
            }
        }

        numGenerations++;
        attempts++;
    }//end while

    //reestablish the original values
    Parameters.setMinCoveredExamplesRule(auxMinCov);
    Parameters.setMaxCoveredExamplesRule(auxMaxCov);
    Parameters.setSizeTournament(auxTournament);
}//end constructor

public ArrayList<Individual> getCurrentPopulation() {
    return currentPopulation;
}

public Individual getBestIndividual() {
    return bestIndividual;
}

}
