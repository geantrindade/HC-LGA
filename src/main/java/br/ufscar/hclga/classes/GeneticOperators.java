package br.ufscar.hclga.classes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class GeneticOperators {

static boolean multiObj;

/* ===========================================================
     * Apply the elitism operator
     * =========================================================== */
public static ArrayList<Individual> elitism(Evolution e) {
    ArrayList<Individual> elite = new ArrayList<Individual>();

    for (int i = 0; i < Parameters.getElitismNumber(); i++) {
        elite.add(e.getCurrentPopulation().get(i));
    }

    return elite;
}

/* ===========================================================
     * Apply the uniform crossover operator considering distance
     * between the parents
     * =========================================================== */
public static ArrayList<Individual> uniformCrossoverDistance(Evolution e) {
    ArrayList<Individual> parents = new ArrayList<Individual>();
    ArrayList<Individual> children = new ArrayList<Individual>();

    int numParents = e.getCurrentPopulation().size() - Parameters.getElitismNumber();
    int numChildren = numParents;

    if (numParents % 2 != 0) {
        numParents++;
    }

    for (int i = 0; i < numParents; i++) {
        parents.add(tournamentSelection1Parent(e));
    }

    while (parents.size() > 0) {
        Individual parent1 = parents.get(0);
        parents.remove(0);

        Individual parent2 = getSecondParent(parent1, parents);

        //Do uniform crossover according to the crossover rate
        Random generator = new Random();

        if (generator.nextDouble() < Parameters.getCrossoverRate()) {
            double[] rule1 = parent1.getRule().clone();
            double[] rule2 = parent2.getRule().clone();

            int numActiveTermsRule1 = 0;
            int numActiveTermsRule2 = 0;

            ArrayList<Integer> posActiveTermsRule1 = new ArrayList<Integer>();
            ArrayList<Integer> posActiveTermsRule2 = new ArrayList<Integer>();
            ArrayList<Integer> activeTerms = new ArrayList<Integer>();

            for (int pos = 0; pos < parent1.getPosActiveTerms().size(); pos++) {
                activeTerms.add(parent1.getPosActiveTerms().get(pos));
            }

            activeTerms = addAllWithoutRepetitions(activeTerms, parent2.getPosActiveTerms());
            int j;

            for (int k = 0; k < activeTerms.size(); k++) {
                j = activeTerms.get(k);

                if (generator.nextDouble() < 0.5) {
                    //exchange the antecedents
                    rule1[j] = parent2.getRule()[j];
                    numActiveTermsRule1 += rule1[j];

                    if (rule1[j] == 1) {
                        posActiveTermsRule1.add(j);
                    }

                    rule1[j + 1] = parent2.getRule()[j + 1];
                    rule1[j + 2] = parent2.getRule()[j + 2];
                    rule1[j + 3] = parent2.getRule()[j + 3];

                    rule2[j] = parent1.getRule()[j];
                    numActiveTermsRule2 += rule2[j];

                    if (rule2[j] == 1) {
                        posActiveTermsRule2.add(j);
                    }

                    rule2[j + 1] = parent1.getRule()[j + 1];
                    rule2[j + 2] = parent1.getRule()[j + 2];
                    rule2[j + 3] = parent1.getRule()[j + 3];

                } else {
                    //copy the antecedents
                    rule1[j] = parent1.getRule()[j];
                    numActiveTermsRule1 += rule1[j];

                    if (rule1[j] == 1) {
                        posActiveTermsRule1.add(j);
                    }

                    rule1[j + 1] = parent1.getRule()[j + 1];
                    rule1[j + 2] = parent1.getRule()[j + 2];
                    rule1[j + 3] = parent1.getRule()[j + 3];

                    rule2[j] = parent2.getRule()[j];
                    numActiveTermsRule2 += rule2[j];

                    if (rule2[j] == 1) {
                        posActiveTermsRule2.add(j);
                    }

                    rule2[j + 1] = parent2.getRule()[j + 1];
                    rule2[j + 2] = parent2.getRule()[j + 2];
                    rule2[j + 3] = parent2.getRule()[j + 3];
                }
            }

            //I will not allow rules with no terms.
            if (numActiveTermsRule1 == 0 || numActiveTermsRule2 == 0) {
                parents.add(parent1);
                parents.add(parent2);

            } else {
                Individual child1 = new Individual(rule1, posActiveTermsRule1, multiObj);
                Individual child2 = new Individual(rule2, posActiveTermsRule2, multiObj);

                if (children.size() < numChildren) {
                    children.add(child1);
                }

                if (children.size() < numChildren) {
                    children.add(child2);
                }
            }

        } else {
            Individual child1 = new Individual(parent1.getRule(), parent1.getPosActiveTerms(), multiObj);
            Individual child2 = new Individual(parent2.getRule(), parent2.getPosActiveTerms(), multiObj);

            if (children.size() < numChildren) {
                children.add(child1);
            }

            if (children.size() < numChildren) {
                children.add(child2);
            }

        }
    }
    return children;
}

/* ===========================================================
     * Apply the uniform crossover operator
     * =========================================================== */
public static ArrayList<Individual> uniformCrossoverRandom(Evolution e) {
    ArrayList<Individual> children = new ArrayList<Individual>();
    int numChildren = e.getCurrentPopulation().size() - Parameters.getElitismNumber();

    while (children.size() < numChildren) {
        ArrayList<Individual> parents = new ArrayList<Individual>();
        parents = tournamentSelection2Parents(e);

        Individual parent1 = parents.get(0);
        Individual parent2 = parents.get(1);

        //Do uniform crossover according to the croosver rate
        Random generator = new Random();

        if (generator.nextDouble() < Parameters.getCrossoverRate()) {
            int numActiveTermsRule1 = 0;
            int numActiveTermsRule2 = 0;

            double[] rule1 = parent1.getRule().clone();
            double[] rule2 = parent2.getRule().clone();

            ArrayList<Integer> posActiveTermsRule1 = new ArrayList<Integer>();
            ArrayList<Integer> posActiveTermsRule2 = new ArrayList<Integer>();
            ArrayList<Integer> activeTerms = new ArrayList<Integer>();

            for (int pos = 0; pos < parent1.getPosActiveTerms().size(); pos++) {
                activeTerms.add(parent1.getPosActiveTerms().get(pos));
            }

            activeTerms = addAllWithoutRepetitions(activeTerms, parent2.getPosActiveTerms());
            int i;

            for (int k = 0; k < activeTerms.size(); k++) {
                i = activeTerms.get(k);

                if (generator.nextDouble() < 0.5) {
                    //exchange the antecedents
                    rule1[i] = parent2.getRule()[i];
                    numActiveTermsRule1 += rule1[i];

                    if (rule1[i] == 1) {
                        posActiveTermsRule1.add(i);
                    }

                    rule1[i + 1] = parent2.getRule()[i + 1];
                    rule1[i + 2] = parent2.getRule()[i + 2];
                    rule1[i + 3] = parent2.getRule()[i + 3];

                    rule2[i] = parent1.getRule()[i];
                    numActiveTermsRule2 += rule2[i];

                    if (rule2[i] == 1) {
                        posActiveTermsRule2.add(i);
                    }

                    rule2[i + 1] = parent1.getRule()[i + 1];
                    rule2[i + 2] = parent1.getRule()[i + 2];
                    rule2[i + 3] = parent1.getRule()[i + 3];

                } else {
                    //copy the antecedents
                    rule1[i] = parent1.getRule()[i];
                    numActiveTermsRule1 += rule1[i];

                    if (rule1[i] == 1) {
                        posActiveTermsRule1.add(i);
                    }

                    rule1[i + 1] = parent1.getRule()[i + 1];
                    rule1[i + 2] = parent1.getRule()[i + 2];
                    rule1[i + 3] = parent1.getRule()[i + 3];

                    rule2[i] = parent2.getRule()[i];
                    numActiveTermsRule2 += rule2[i];

                    if (rule2[i] == 1) {
                        posActiveTermsRule2.add(i);
                    }

                    rule2[i + 1] = parent2.getRule()[i + 1];
                    rule2[i + 2] = parent2.getRule()[i + 2];
                    rule2[i + 3] = parent2.getRule()[i + 3];
                }
            }

            Individual child1 = new Individual(rule1, posActiveTermsRule1, multiObj);
            Individual child2 = new Individual(rule2, posActiveTermsRule2, multiObj);

            if (numActiveTermsRule1 > 0 && numActiveTermsRule2 > 0) {
                if (children.size() < numChildren) {
                    children.add(child1);
                }

                if (children.size() < numChildren) {
                    children.add(child2);
                }
            }

        } else {

            Individual child1 = new Individual(parent1.getRule(), parent1.getPosActiveTerms(), multiObj);
            Individual child2 = new Individual(parent2.getRule(), parent2.getPosActiveTerms(), multiObj);

            if (children.size() < numChildren) {
                children.add(child1);
            }

            if (children.size() < numChildren) {
                children.add(child2);
            }
        }

    }

    return children;
}

/* ===========================================================
     * Given a rule, returns the position of the active terms
     * =========================================================== */
public static ArrayList<Integer> getPosActiveTerms(double[] rule1, double[] rule2) {
    ArrayList<Integer> posActiveTerms = new ArrayList<Integer>();

    for (int i = 0; i < rule1.length; i += 4) {
        if (rule1[i] == 1 && posActiveTerms.contains(i) == false) {
            posActiveTerms.add(i);
        }

        if (rule2[i] == 1 && posActiveTerms.contains(i) == false) {
            posActiveTerms.add(i);
        }
    }

    return posActiveTerms;
}

public static ArrayList<Integer> getPosActiveTerms(double[] rule) {
    ArrayList<Integer> posActiveTerms = new ArrayList<Integer>();

    for (int i = 0; i < rule.length; i += 4) {
        if (rule[i] == 1) {
            posActiveTerms.add(i);
        }
    }

    return posActiveTerms;
}

/* ===========================================================
     * Select one parent via tournament selection
     * =========================================================== */
private static Individual tournamentSelection1Parent(Evolution e) {
    int sizeTournament = Parameters.getSizeTournament();

    Random generator = new Random();
    ArrayList<Individual> candidates = new ArrayList<Individual>();

    //Get the candidates
    int i = 0;

    while (i < sizeTournament) {
        int index = generator.nextInt(e.getCurrentPopulation().size());

        if (candidates.contains(e.getCurrentPopulation().get(index)) == false) {
            candidates.add(e.getCurrentPopulation().get(index));
            i++;
        }
    }

    Collections.sort(candidates);
    //Return the best candidate
    return candidates.get(0);
}

/* ===========================================================
     * Select two parents via tournament selection
     * =========================================================== */
private static ArrayList<Individual> tournamentSelection2Parents(Evolution e) {
    ArrayList<Individual> candidates = new ArrayList<Individual>();
    int sizeTournament = Parameters.getSizeTournament();

    Random generator = new Random();
    ArrayList<Individual> parents = new ArrayList<Individual>();

    //Get ehe first parent
    candidates.clear();
    int i = 0;

    while (i < sizeTournament) {
        int index = generator.nextInt(e.getCurrentPopulation().size());

        if (candidates.contains(e.getCurrentPopulation().get(index)) == false) {
            candidates.add(e.getCurrentPopulation().get(index));
            i++;
        }
    }

    Collections.sort(candidates);
    parents.add(candidates.get(0));

    //Get the second parent
    candidates.clear();
    i = 0;

    while (i < sizeTournament) {
        int index = generator.nextInt(e.getCurrentPopulation().size());

        if (candidates.contains(e.getCurrentPopulation().get(index)) == false
                && parents.contains(e.getCurrentPopulation().get(index)) == false) {
            candidates.add(e.getCurrentPopulation().get(index));
            i++;
        }
    }

    Collections.sort(candidates);
    parents.add(candidates.get(0));

    return parents;
}

/* ===========================================================
     * Get the second parent for crossover
     * =========================================================== */
private static Individual getSecondParent(Individual parent1, ArrayList<Individual> parents) {
    Individual parent2;
    ArrayList<Distances> distances = new ArrayList<Distances>();

    for (int i = 0; i < parents.size(); i++) {
        Individual candidate = parents.get(i);
        ArrayList<Integer> posActiveTerms1 = parent1.getPosActiveTerms();
        ArrayList<Integer> posActiveTerms2 = candidate.getPosActiveTerms();

        if (posActiveTerms1.equals(posActiveTerms2) == false) {
            distances.add(new Distances(getEuclideanDistance(parent1, candidate), i));
        }
    }

    Collections.sort(distances);

    if (distances.isEmpty()) {
        parent2 = parents.get(0);
        int indexParent2 = 0;
        parents.remove(indexParent2);

    } else {
        Distances shortestDistance = distances.get(0);
        int indexParent2 = shortestDistance.getPosition();
        parent2 = parents.get(indexParent2);
        parents.remove(indexParent2);
    }

    return parent2;
}

/* ===========================================================
     * Get the Euclidean distance between two mean class label
     * vectors of two individuals
     * =========================================================== */
private static double getEuclideanDistance(Individual rule1, Individual rule2) {
    double distance = 0;

    double[] meanClassLabel1 = rule1.getMeanClassLabelVectorCovered();
    double[] meanClassLabel2 = rule2.getMeanClassLabelVectorCovered();
    double[] weights = Classes.getWeightingScheme();

    double sum = 0;

    for (int i = 0; i < meanClassLabel1.length; i++) {
        sum += weights[i] * Math.pow((meanClassLabel1[i] - meanClassLabel2[i]), 2);
    }

    distance = Math.sqrt(sum);

    return distance;
}

/* ===========================================================
     * Tries to search for a rule with good fitness and with a
     * maximum coverage
     * =========================================================== */
public static ArrayList<Individual> localSearchOperatorMaxFitness(ArrayList<Individual> individuals) {
    Random generator = new Random();

    for (int i = 0; i < individuals.size(); i++) {
        Individual individual = individuals.get(i);
        double[] previousRule = new double[individual.getRule().length];
        ArrayList<Integer> previousPosActiveTerms = new ArrayList<Integer>();
        double previousFitness = 0;
        int first = 0;

        while ((individual.getFitness() >= previousFitness
                && individual.getNumberCoveredExamples() < Parameters.getMaxCoveredExamplesRule()) || first == 0) {

            first++;

            //Save the rule before changing
            previousRule = individual.getRule().clone();
            previousFitness = individual.getFitness();
            previousPosActiveTerms = individual.getPosActiveTerms();

            if (individual.getPosActiveTerms().size() > 1) {
                //Let's remove one term
                int pos = generator.nextInt(individual.getPosActiveTerms().size());
                int posToMutate = individual.getPosActiveTerms().get(pos);
                individual.getRule()[posToMutate] = 0.0;
                individual.getPosActiveTerms().remove(pos);

            } else {
                //Randomly substitute the actual term by another
                individual.getRule()[individual.getPosActiveTerms().get(0)] = 0.0;
                individual.getPosActiveTerms().remove(0);

                int posTerm = generator.nextInt(individual.getRule().length);

                while ((posTerm % 4) != 0) {
                    posTerm = generator.nextInt(individual.getRule().length);
                }

                individual.getRule()[posTerm] = 1.0;
                individual.getPosActiveTerms().add(posTerm);
            }
            //Substitute the individual by the mutated
            individual = new Individual(individual.getRule(), individual.getPosActiveTerms(), multiObj);
            individuals.set(i, individual);
        }
        //Return the previous individual
        individual = new Individual(previousRule, previousPosActiveTerms, multiObj);
        individuals.set(i, individual);
    }

    return individuals;
}

/* ===========================================================
     * Try to garantee that a rule cover at least a minimum number of
     * examples, and also garantees that a rule don't cover more
     * than a maximum number of examples
     * =========================================================== */
public static ArrayList<Individual> localSearchOperatorMinMax(ArrayList<Individual> individuals) {
    Random generator = new Random();

    for (int i = 0; i < individuals.size(); i++) {
        Individual individual = individuals.get(i);
        int maxAttempts = 0;
        int convergence = 0;

        while (convergence == 0 && maxAttempts < 100) {
            //Will try 100 times
            maxAttempts++;

            if (individual.getNumberCoveredExamples() < Parameters.getMinCoveredExamplesRule()) {
                //The rule is too specific. Let's remove one term
                if (individual.getPosActiveTerms().size() > 1) {
                    int pos = generator.nextInt(individual.getPosActiveTerms().size());
                    int posToMutate = individual.getPosActiveTerms().get(pos);
                    individual.getRule()[posToMutate] = 0.0;
                    individual.getPosActiveTerms().remove(pos);

                } else if (individual.getPosActiveTerms().size() == 1) {
                    //Randomly substitute the actual term by another and generalize this term
                    individual.getRule()[individual.getPosActiveTerms().get(0)] = 0.0;
                    individual.getPosActiveTerms().remove(0);

                    int posTerm = generator.nextInt(individual.getRule().length);
                    while ((posTerm % 4) != 0) {
                        posTerm = generator.nextInt(individual.getRule().length);
                    }

                    individual.getRule()[posTerm] = 1.0;
                    individual.getPosActiveTerms().add(posTerm);

                    int posAttribute = posTerm / 4;
                    if (Datasets.getInfoAttributes().get(posAttribute) == 1) {
                        generalizationNumeric(individual.getRule(), posTerm, 1);
                    } else {
                        int posCatAttribute = Datasets.getCategoricAttributes().indexOf(posAttribute);
                        generalizationCategorical(individual.getRule(), posTerm, posCatAttribute);
                    }
                } else {
                    //Empty rule!!!!
                    //Add one term
                    int posTerm = generator.nextInt(individual.getRule().length);

                    while ((posTerm % 4) != 0) {
                        posTerm = generator.nextInt(individual.getRule().length);
                    }

                    individual.getRule()[posTerm] = 1.0;
                    individual.getPosActiveTerms().add(posTerm);
                }
                //Substitute the individual by the mutated
                individual = new Individual(individual.getRule(), individual.getPosActiveTerms(), multiObj);
                individuals.set(i, individual);

            } else if (individual.getNumberCoveredExamples() > Parameters.getMaxCoveredExamplesRule()) {
                //The rule is too general. Let's add one term
                if (individual.getPosActiveTerms().size() == Datasets.getInfoAttributes().size()) {
                    //Results.printRule(individual.getRule());
                    //Wow, general rule that uses all atributes. Let's give up
                    //and let evolution do its job
                    convergence = 1;

                } else {
                    int posToMutate = generator.nextInt(Datasets.getInfoAttributes().size()) * 4;

                    while (individual.getPosActiveTerms().contains(posToMutate) == true) {
                        posToMutate = generator.nextInt(Datasets.getInfoAttributes().size()) * 4;
                    }

                    individual.getRule()[posToMutate] = 1.0;
                    individual.getPosActiveTerms().add(posToMutate);

                    //Substitute the individual by the mutated
                    individual = new Individual(individual.getRule(), individual.getPosActiveTerms(), multiObj);
                    individuals.set(i, individual);
                }

            } else {
                convergence = 1;
            }
        }
    }

    return individuals;
}

/* ===========================================================
     * Eliminates terms while the fitness function does not decrease
     * =========================================================== */
public static ArrayList<Individual> localSearchOperatorFitness(ArrayList<Individual> individuals) {
    Individual individual;
    double previousFitness;

    ArrayList<Integer> activeTerms = new ArrayList<Integer>();

    for (int i = 0; i < individuals.size(); i++) {
        individual = individuals.get(i);
        previousFitness = individual.getFitness();

        for (int pos = 0; pos < individual.getPosActiveTerms().size(); pos++) {
            activeTerms.add(individual.getPosActiveTerms().get(pos));
        }

        for (int j = 0; j < activeTerms.size(); j++) {
            individual.getRule()[activeTerms.get(j)] = 0.0;
            individual.getPosActiveTerms().remove(activeTerms.get(j));
            individual = new Individual(individual.getRule(), individual.getPosActiveTerms(), multiObj);

            if (individual.getFitness() < previousFitness) {
                individual.getRule()[activeTerms.get(j)] = 1.0;
                individual.getPosActiveTerms().add(activeTerms.get(j));
                individual = new Individual(individual.getRule(), individual.getPosActiveTerms(), multiObj);

            } else {
                previousFitness = individual.getFitness();
            }
        }

        individuals.set(i, individual);
    }

    return individuals;
}

/* ===========================================================
     * Apply the mutation operator
     * =========================================================== */
public static ArrayList<Individual> mutation(ArrayList<Individual> children) {
    int numMutationIndividuals = (int) (Parameters.getMutationRate() * children.size());
    ArrayList<Integer> indexChildrenMutation = new ArrayList<Integer>();

    //Get the individuals that will suffer mutation
    Random generator = new Random();
    for (int i = 0; i < numMutationIndividuals; i++) {
        int index = generator.nextInt(children.size());
        if (indexChildrenMutation.contains(index) == false) {
            indexChildrenMutation.add(index);
        }
    }

    //Apply the mutations
    for (int i = 0; i < indexChildrenMutation.size(); i++) {
        int indexIndividual = indexChildrenMutation.get(i);
        ArrayList<Integer> posActiveTerms = new ArrayList<Integer>();

        //Flag mutation?
        if (generator.nextDouble() <= 0.5) {
            //Change the flags with Parameters.getProbabilityUseClausule() probability
            // Will not allow the mutation to generate empty rules
            double[] rule = new double[children.get(indexIndividual).getRule().length];
            System.arraycopy(children.get(indexIndividual).getRule(), 0, rule, 0, rule.length);
            int numActives = 0;

            while (numActives == 0) {
                posActiveTerms = new ArrayList<Integer>();

                for (int posTerm = 0; posTerm < children.get(indexIndividual).getRule().length; posTerm += 4) {
                    numActives += flagMutation(rule, posTerm, Parameters.getProbabilityUseClausule(), posActiveTerms);
                    //numActives += flagMutation(rule, posTerm, 0.5, posActiveTerms);
                }
            }

            Individual individual = new Individual(rule, posActiveTerms, multiObj);
            children.set(indexIndividual, individual);
        } //No flag mutation?
        else {
            //Generalize or restrict the rules with 0.5 probability
            for (int j = 0; j < children.get(indexIndividual).getPosActiveTerms().size(); j++) {
                int posTerm = children.get(indexIndividual).getPosActiveTerms().get(j);
                //Ramdomly choose an mutation operation depending on the attribute
                //value, if categorical or numeric
                int posAttribute = posTerm / 4;

                if (Datasets.getInfoAttributes().get(posAttribute) == 1) {
                    //Numeric attribute, so choose between generalization and restriction
                    //0 -> Generalization
                    //1 -> Restriction
                    int operation = generator.nextInt(2);

                    switch (operation) {
                    case 0:
                        generalizationNumeric(children.get(indexIndividual).getRule(), posTerm, 1);
                        break;

                    case 1:
                        restrictionNumeric(children.get(indexIndividual).getRule(), posTerm, 1);
                        break;
                    }
                } else {
                    //Categoric attribute, so apply generalization or restriction
                    //depending on the values of the rule term
                    int posCatAttribute = Datasets.getCategoricAttributes().indexOf(posAttribute);
                    generalizationRestrictionCategorical(children.get(indexIndividual).getRule(), posTerm, posCatAttribute, 1);
                }
            }
            //Substitute the individual by the mutated
            Individual individual = new Individual(children.get(indexIndividual).getRule(), children.get(indexIndividual).getPosActiveTerms(), multiObj);
            children.set(indexIndividual, individual);
        }
    }

    return children;
}

/* ===========================================================
     * Generalization for categoric attribute
     * =========================================================== */
private static void generalizationCategorical(double[] rule, int posTerm, int posAttribute) {
    Random generator = new Random();

    //Get the operator of the term in the rule
    double ruleOperator = rule[posTerm + 1];

    if (ruleOperator == 0) { //attributeValue = infLim
        //Term cannot be restricted, so we will randomly substitute the "=" operator by "in" or "!="
        if (generator.nextDouble() < 0.5) {
            //Substitute the "=" operator by the "!=" operator
            rule[posTerm + 1] = 1.0;

        } else {
            //Substitute the "=" operator by the "in" operator, which adds
            //more categorical values to the condition, generalizing the rule
            rule[posTerm + 1] = 2.0;

            //Get the categorical attribute value
            double infLim = rule[posTerm + 2];
            String attributeValue = Datasets.getCategoricMapping().get(posAttribute).get((int) infLim)[0];

            //Will store the indexes of the CategoricMapping that can substitute
            //the actual value
            ArrayList<Double> possibleSubstitutions = new ArrayList<Double>();

            int numPossibleValues = Datasets.getCategoricMapping().get(posAttribute).size();

            for (int i = 0; i < numPossibleValues; i++) {
                if (Datasets.getCategoricMapping().get(posAttribute).get(i).length > 1) {
                    for (int j = 0; j < Datasets.getCategoricMapping().get(posAttribute).get(i).length; j++) {
                        if (Datasets.getCategoricMapping().get(posAttribute).get(i)[j].equals(attributeValue) == true) {
                            possibleSubstitutions.add((double) i);
                            break;
                        }
                    }
                }
            }
            //Randomly choose a substitution
            int newIndexValue = generator.nextInt(possibleSubstitutions.size());
            rule[posTerm + 2] = possibleSubstitutions.get(newIndexValue);
        }
    } else if (ruleOperator == 2) { //attributeValue in {values in infLim}
        //Will randomly substitute the set of values of the term by a 
        //larger set. 
        //Get the size of the actual set of value in the term
        double infLim = rule[posTerm + 2];
        String[] setActualValues = Datasets.getCategoricMapping().get(posAttribute).get((int) infLim);

        //Will store the indexes of the CategoricMapping that can substitute
        //the actual value
        ArrayList<Double> possibleSubstitutions = new ArrayList<Double>();

        int numPossibleValues = Datasets.getCategoricMapping().get(posAttribute).size();

        for (int i = 0; i < numPossibleValues; i++) {
            if (Arrays.equals(Datasets.getCategoricMapping().get(posAttribute).get(i), setActualValues) == false
                    && Datasets.getCategoricMapping().get(posAttribute).get(i).length > setActualValues.length) {
                possibleSubstitutions.add((double) i);
            }
        }
        //If the actual set of values is already the bigger, cannot generalize. :-(
        if (possibleSubstitutions.isEmpty() == false) {
            //Randomly choose a substitution
            int newIndexValue = generator.nextInt(possibleSubstitutions.size());
            rule[posTerm + 2] = possibleSubstitutions.get(newIndexValue);
        }
    }
}

/* ===========================================================
     * Generalization/Restriction mutation for categoric attribute
     * =========================================================== */
private static void generalizationRestrictionCategorical(double[] rule, int posTerm, int posAttribute, int probability) {
    //Get the operator of the term in the rule
    double ruleOperator = rule[posTerm + 1];

    Random generator = new Random();
    double prob = generator.nextDouble();

    if (prob <= probability) {
        if (ruleOperator == 0) { //attributeValue = infLim
            //Term cannot be restricted, so we will randomly substitute the "=" operator by "in" or "!="
            if (generator.nextDouble() < 0.5) {
                //Substitute the "=" operator by the "!=" operator
                rule[posTerm + 1] = 1.0;

            } else {
                //Substitute the "=" operator by the "in" operator, which adds
                //more categorical values to the condition, generalizing the rule
                rule[posTerm + 1] = 2.0;

                //Get the categorical attribute value
                double infLim = rule[posTerm + 2];
                String attributeValue = Datasets.getCategoricMapping().get(posAttribute).get((int) infLim)[0];

                //Will store the indexes of the CategoricMapping that can substitute
                //the actual value
                ArrayList<Double> possibleSubstitutions = new ArrayList<Double>();

                int numPossibleValues = Datasets.getCategoricMapping().get(posAttribute).size();

                for (int i = 0; i < numPossibleValues; i++) {
                    if (Datasets.getCategoricMapping().get(posAttribute).get(i).length > 1) {
                        for (int j = 0; j < Datasets.getCategoricMapping().get(posAttribute).get(i).length; j++) {
                            if (Datasets.getCategoricMapping().get(posAttribute).get(i)[j].equals(attributeValue) == true) {
                                possibleSubstitutions.add((double) i);
                                break;
                            }
                        }
                    }
                }
                //Randomly choose a substitution
                int newIndexValue = generator.nextInt(possibleSubstitutions.size());
                rule[posTerm + 2] = possibleSubstitutions.get(newIndexValue);
            }
        } else if (ruleOperator == 1) { //attributeValue != infLim
            //Term can be generalized or restricted depending on the attribute values
            //So we will randomly substitute the "!=" operator by "in" or "="
            if (generator.nextDouble() < 0.5) {
                //Substitute the "!=" operator by the "=" operator
                rule[posTerm + 1] = 0.0;

            } else {
                //Substitute the "!=" operator by the "in" operator
                //Here the rule will be randomly restricted or generalized
                //We will randomly substitute the value in the term by a set of values
                rule[posTerm + 1] = 2.0;

                //Will store the indexes of the CategoricMapping that can substitute
                //the actual value
                ArrayList<Double> possibleSubstitutions = new ArrayList<Double>();

                int numPossibleValues = Datasets.getCategoricMapping().get(posAttribute).size();

                for (int i = 0; i < numPossibleValues; i++) {
                    if (Datasets.getCategoricMapping().get(posAttribute).get(i).length > 1) {
                        possibleSubstitutions.add((double) i);
                    }
                }

                //Randomly choose a substitution
                int newIndexValue = generator.nextInt(possibleSubstitutions.size());
                rule[posTerm + 2] = possibleSubstitutions.get(newIndexValue);
            }
        } else if (ruleOperator == 2) { //attributeValue in {values in infLim}
            //Will randomly substitute the set of values of the term by a smaller
            //or larger set. If the smaller set contains just one element, will substitute
            //the "in" operator by the "=" operator
            //Get the size of the actual set of value in the term
            double infLim = rule[posTerm + 2];
            String[] setActualValues = Datasets.getCategoricMapping().get(posAttribute).get((int) infLim);

            //Will store the indexes of the CategoricMapping that can substitute
            //the actual value
            ArrayList<Double> possibleSubstitutions = new ArrayList<Double>();

            int numPossibleValues = Datasets.getCategoricMapping().get(posAttribute).size();

            for (int i = 0; i < numPossibleValues; i++) {
                if (Arrays.equals(Datasets.getCategoricMapping().get(posAttribute).get(i), setActualValues) == false) {
                    possibleSubstitutions.add((double) i);
                }
            }
            //Randomly choose a substitution
            int newIndexValue = generator.nextInt(possibleSubstitutions.size());
            rule[posTerm + 2] = possibleSubstitutions.get(newIndexValue);
            infLim = rule[posTerm + 2];
            //If the substitution set of values contains just one element,
            //substitute the "in" operator by the "=" operator
            if (Datasets.getCategoricMapping().get(posAttribute).get((int) infLim).length == 1) {
                rule[posTerm + 1] = 0.0;
            }
        }
    }
}

/* ===========================================================
     * Generalization mutation for numeric attribute
     * =========================================================== */
private static void generalizationNumeric(double[] rule, int posTerm, double probability) {

    //Get the operator of the term in the rule
    double ruleOperator = rule[posTerm + 1];

    Random generator = new Random();
    double prob = generator.nextDouble();

    if (prob <= probability) {

        double num = generator.nextDouble();

        if (ruleOperator == 0) { //Attribute <= supLim
            double factor = rule[posTerm + 3] * num;

            if (rule[posTerm + 3] == 0) {
                rule[posTerm + 3] += num;
            } else if (rule[posTerm + 3] > 0) {
                rule[posTerm + 3] += factor;
            } else {
                rule[posTerm + 3] -= factor;
            }

        } else if (ruleOperator == 1) { //Attribute >= infLim
            double factor = rule[posTerm + 2] * num;

            if (rule[posTerm + 2] == 0) {
                rule[posTerm + 2] -= num;
            } else if (rule[posTerm + 2] > 0) {
                rule[posTerm + 2] -= factor;
            } else {
                rule[posTerm + 2] += factor;
            }

        } else if (ruleOperator == 2) {//infLim <= Attribute <= supLim
            double factorInf = rule[posTerm + 2] * num;
            double factorSup = rule[posTerm + 3] * num;

            if (rule[posTerm + 2] == 0) {
                rule[posTerm + 2] -= num;
            } else if (rule[posTerm + 2] > 0) {
                rule[posTerm + 2] -= factorInf;
            } else {
                rule[posTerm + 2] += factorInf;
            }

            if (rule[posTerm + 3] == 0) {
                rule[posTerm + 3] += num;
            } else if (rule[posTerm + 3] > 0) {
                rule[posTerm + 3] += factorSup;
            } else {
                rule[posTerm + 3] -= factorSup;
            }
        }
    }
}

/* ===========================================================
     * Restriction mutation for numeric attribute
     * =========================================================== */
private static void restrictionNumeric(double[] rule, int posTerm, double probability) {
    //Get the operator of the term in the rule
    double ruleOperator = rule[posTerm + 1];

    Random generator = new Random();
    double prob = generator.nextDouble();

    if (prob <= probability) {
        double num = generator.nextDouble();

        if (ruleOperator == 0) { //Attribute <= supLim
            double factor = rule[posTerm + 3] * num;

            if (rule[posTerm + 3] == 0) {
                rule[posTerm + 3] -= num;
            } else if (rule[posTerm + 3] > 0) {
                rule[posTerm + 3] -= factor;
            } else {
                rule[posTerm + 3] += factor;
            }

        } else if (ruleOperator == 1) { //Attribute >= infLim
            double factor = rule[posTerm + 2] * num;

            if (rule[posTerm + 2] == 0) {
                rule[posTerm + 2] += num;
            } else if (rule[posTerm + 2] > 0) {
                rule[posTerm + 2] += factor;
            } else {
                rule[posTerm + 2] -= factor;
            }

        } else if (ruleOperator == 2) {//infLim <= Attribute <= supLim
            double factorInf = rule[posTerm + 2] * num;
            double factorSup = rule[posTerm + 3] * num;

            if (rule[posTerm + 2] == 0) {
                rule[posTerm + 2] += num;
            } else if (rule[posTerm + 2] > 0) {
                rule[posTerm + 2] += factorInf;
            } else {
                rule[posTerm + 2] -= factorInf;
            }

            if (rule[posTerm + 3] == 0) {
                rule[posTerm + 3] -= num;
            } else if (rule[posTerm + 3] > 0) {
                rule[posTerm + 3] -= factorSup;
            } else {
                rule[posTerm + 3] += factorSup;
            }
        }
    }
}

/* ===========================================================
     * Flag mutation
     * =========================================================== */
public static int flagMutation(Individual individual, int pos, double probability, ArrayList<Integer> posActiveTerms) {
    Random generator = new Random();
    double num = generator.nextDouble();
    int active = 0;

    if (num <= probability) {
        if (individual.getRule()[pos] == 0) {
            individual.getRule()[pos] = 1.0;
            active = 1;
            posActiveTerms.add(pos);
        } else {
            individual.getRule()[pos] = 0.0;
        }

    } else if (individual.getRule()[pos] == 1) {
        active = 1;
        posActiveTerms.add(pos);
    }

    return active;
}

public static int flagMutation(double rule[], int pos, double probability, ArrayList<Integer> posActiveTerms) {
    Random generator = new Random();
    double num = generator.nextDouble();
    int active = 0;

    if (num <= probability) {
        if (rule[pos] == 0) {
            rule[pos] = 1.0;
            active = 1;
            posActiveTerms.add(pos);

        } else {
            rule[pos] = 0.0;
        }

    } else if (rule[pos] == 1) {
        active = 1;
        posActiveTerms.add(pos);
    }

    return active;
}

/* ===========================================================
     * Add one array list into another without repetions
     * =========================================================== */
public static ArrayList<Integer> addAllWithoutRepetitions(ArrayList<Integer> terms, ArrayList<Integer> termsToAdd) {
    for (int i = 0; i < termsToAdd.size(); i++) {
        if (terms.contains(termsToAdd.get(i)) == false) {
            terms.add(termsToAdd.get(i));
        }
    }

    return terms;
}

}

/* ===========================================================
 * Class to order the Euclidean distances
 * =========================================================== */
class Distances implements Comparable {

private double distance;
private int position;

public Distances(double dist, int pos) {
    distance = dist;
    position = pos;
}

/* ===========================================================
     * Sorts the distances
     * =========================================================== */
public int compareTo(Object o) {
    if (this.getDistance() < ((Distances) o).getDistance()) {
        return -1;
    } else if (this.getDistance() > ((Distances) o).getDistance()) {
        return 1;
    } else {
        return 0;
    }
}

public double getDistance() {
    return distance;
}

public int getPosition() {
    return position;
}
}
