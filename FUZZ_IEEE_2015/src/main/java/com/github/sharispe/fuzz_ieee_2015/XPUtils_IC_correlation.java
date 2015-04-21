/*
 *  Copyright or © or Copr. Ecole des Mines d'Alès (2012-2014) 
 *  
 *  This software is a computer program whose purpose is to provide 
 *  several functionalities for the processing of semantic data 
 *  sources such as ontologies or text corpora.
 *  
 *  This software is governed by the CeCILL  license under French law and
 *  abiding by the rules of distribution of free software.  You can  use, 
 *  modify and/ or redistribute the software under the terms of the CeCILL
 *  license as circulated by CEA, CNRS and INRIA at the following URL
 *  "http://www.cecill.info". 
 * 
 *  As a counterpart to the access to the source code and  rights to copy,
 *  modify and redistribute granted by the license, users are provided only
 *  with a limited warranty  and the software's author,  the holder of the
 *  economic rights,  and the successive licensors  have only  limited
 *  liability. 

 *  In this respect, the user's attention is drawn to the risks associated
 *  with loading,  using,  modifying and/or developing or reproducing the
 *  software by the user in light of its specific status of free software,
 *  that may mean  that it is complicated to manipulate,  and  that  also
 *  therefore means  that it is reserved for developers  and  experienced
 *  professionals having in-depth computer knowledge. Users are therefore
 *  encouraged to load and test the software's suitability as regards their
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and,  more generally, to use and operate it in the 
 *  same conditions as regards security. 
 * 
 *  The fact that you are presently reading this means that you have had
 *  knowledge of the CeCILL license and that you accept its terms.
 */
package com.github.sharispe.fuzz_ieee_2015;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openrdf.model.URI;
import slib.graph.model.graph.G;
import slib.graph.model.impl.repo.URIFactoryMemory;
import slib.graph.model.repo.URIFactory;
import slib.sml.sm.core.engine.SM_Engine;
import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Corpus;
import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Topo;
import slib.sml.sm.core.metrics.ic.utils.ICconf;
import slib.sml.sm.core.utils.SMConstants;
import slib.utils.ex.SLIB_Exception;

/**
 *
 * Class used to compute the correlation between several ICs
 *
 * @author Sébastien Harispe <sebastien.harispe@gmail.com>
 */
public class XPUtils_IC_correlation {

    public static void main(String[] args) throws Exception {

        int massFactor = 100;
        boolean useOnlyNouns = true;
        boolean buildConceptRepresentation = false;
        URIFactory factory = URIFactoryMemory.getSingleton();

        G wordNet = DataLoaderUtils.loadWordNet(factory, Conf.dataloc, useOnlyNouns);
        DataLoaderUtils.loadAnnotations(factory, wordNet, Conf.annotations);

        Map<URI, String> labels = DataLoaderUtils.loadURILabels(factory, wordNet, Conf.dataloc, useOnlyNouns);

        FuzzIEEE_XPBrain computer = new FuzzIEEE_XPBrain(factory, wordNet, labels, Conf.file_concept_representation, buildConceptRepresentation, massFactor);

        SM_Engine engine = new SM_Engine(wordNet);

        ICconf icConf_sanchez = new IC_Conf_Topo(SMConstants.FLAG_ICI_SANCHEZ_2011);
        ICconf icConf_seco = new IC_Conf_Topo(SMConstants.FLAG_ICI_SECO_2004);
        ICconf icConf_resnik = new IC_Conf_Corpus(SMConstants.FLAG_IC_ANNOT_RESNIK_1995);

        List<ICconf> confs = new ArrayList();
        confs.add(icConf_sanchez);
        confs.add(icConf_seco);
        confs.add(icConf_resnik);

        List<URI> listConcepts = new ArrayList<>(engine.getClasses());

        System.out.println("Number of concepts: " + listConcepts.size());

        CorrelationComputer corr = new CorrelationComputer();

        int showKmostVariants = 0;
        int nbClasses = 50;

        String report = "mass factor:" + massFactor + ";useOnlyNouns:" + useOnlyNouns + "\n";

        for (CorrelationComputer.COR_METH method : CorrelationComputer.COR_METH.values()) {

            report += "\n----------------------------------------------------\n";
            report += method.toString() + "\n";
            report += "----------------------------------------------------\n";

            for (int k = 0; k < confs.size(); k++) {

                double[] tabk = computeNormalizedICs(engine, confs.get(k), listConcepts);

                for (int l = k + 1; l < confs.size(); l++) {

                    double[] tabl = computeNormalizedICs(engine, confs.get(l), listConcepts);

                    double corrPearson = corr.computeCorrelation(method, tabk, tabl);
                    report += "\n" + confs.get(k).getLabel() + "\t" + confs.get(l).getLabel() + "\t" + corrPearson;

                    report += plotKmostVariants(showKmostVariants, listConcepts, labels, tabk, tabl);

                    System.out.println(confs.get(k).getLabel() + "\t" + confs.get(l).getLabel() + "\t" + corrPearson);

                    computeDiffDistrib(listConcepts, tabk, tabl, nbClasses);
                }
            }

            double[] belief = new double[listConcepts.size()];
            double[] plausibility = new double[listConcepts.size()];
            double[] pignistic = new double[listConcepts.size()];

            // compute belief / plausibility / prob. pignistic for all concepts
            double max_belief = 0, max_plausibility = 0, max_pign = 0;
            for (int k = 0; k < listConcepts.size(); k++) {

                URI concept_k = listConcepts.get(k);
                belief[k] = computer.getIC_Belief(concept_k);
                plausibility[k] = computer.getIC_Plausibility(concept_k);
                pignistic[k] = computer.getIC_PignisticProbabilities(concept_k);

                if (belief[k] > max_belief) {
                    max_belief = belief[k];
                }
                if (plausibility[k] > max_plausibility) {
                    max_plausibility = plausibility[k];
                }
                if (pignistic[k] > max_pign) {
                    max_pign = pignistic[k];
                }
            }

            // norm belief / plausibility / prob. pignistic for all concepts
            for (int k = 0; k < listConcepts.size(); k++) {
                belief[k] = belief[k] / max_belief;
                plausibility[k] = plausibility[k] / max_plausibility;
                pignistic[k] = pignistic[k] / max_pign;
            }

            // belief
            for (int k = 0; k < confs.size(); k++) {
                double[] tabk = computeNormalizedICs(engine, confs.get(k), listConcepts);
                double corrPearson = corr.computeCorrelation(method, tabk, belief);
                report += "\nBelief\t" + confs.get(k).getLabel() + "\t" + corrPearson;

                report += plotKmostVariants(showKmostVariants, listConcepts, labels, tabk, belief);

                System.out.println("\nBelief\t" + confs.get(k).getLabel() + "\t" + corrPearson);
                computeDiffDistrib(listConcepts, tabk, belief, nbClasses);

            }

            // plausibility
            for (int k = 0; k < confs.size(); k++) {
                double[] tabk = computeNormalizedICs(engine, confs.get(k), listConcepts);
                double corrPearson = corr.computeCorrelation(method, tabk, plausibility);
                report += "\nPlausibility\t" + confs.get(k).getLabel() + "\t" + corrPearson;

                report += plotKmostVariants(showKmostVariants, listConcepts, labels, tabk, plausibility);

                System.out.println("\nPlausibility\t" + confs.get(k).getLabel() + "\t" + corrPearson);
                computeDiffDistrib(listConcepts, tabk, plausibility, nbClasses);

            }

            // pignistic
            for (int k = 0; k < confs.size(); k++) {
                double[] tabk = computeNormalizedICs(engine, confs.get(k), listConcepts);
                double corrPearson = corr.computeCorrelation(method, tabk, pignistic);
                report += "\nPignistic\t" + confs.get(k).getLabel() + "\t" + corrPearson;

                report += plotKmostVariants(showKmostVariants, listConcepts, labels, tabk, pignistic);

                System.out.println("\nPignistic\t" + confs.get(k).getLabel() + "\t" + corrPearson);
                computeDiffDistrib(listConcepts, tabk, pignistic, nbClasses);

            }

            double corrPearson_bpl = corr.computeCorrelation(method, belief, plausibility);
            double corrPearson_bpi = corr.computeCorrelation(method, belief, pignistic);
            double corrPearson_plppi = corr.computeCorrelation(method, plausibility, pignistic);
//
            report += "\nBelief\tPlausibility\t" + corrPearson_bpl;
            report += "\nBelief\tPignisti\t" + corrPearson_bpi;
            report += "\nPlausibility\tPignistic\t" + corrPearson_plppi;

        }
        System.out.println(report);

    }

    private static double[] computeNormalizedICs(SM_Engine engine, ICconf icConf, List<URI> listConcepts) throws SLIB_Exception {

        System.out.println("Computing " + icConf.getLabel());

        double max = -1;
        double[] ics = new double[listConcepts.size()];
        for (int i = 0; i < listConcepts.size(); i++) {
            ics[i] = engine.getIC(icConf, listConcepts.get(i));

            if (ics[i] > max) {
                max = ics[i];
            }
        }
        System.out.println("Normalizing");
        for (int i = 0; i < listConcepts.size(); i++) {
            ics[i] = ics[i] / max;
        }

        return ics;
    }

    private static String plotKmostVariants(int nbVariant, List<URI> listConcepts, Map<URI, String> labels, double[] tabk, double[] tabl) {

        if (nbVariant < 1) {
            return "";
        }
        int[] mostVariantICs = getKMostVariant(nbVariant, tabk, tabl);

        String r = "";

        for (int i : mostVariantICs) {

            URI c = listConcepts.get(i);
            double ic_a = tabk[i];
            double ic_b = tabl[i];
            String label = labels.get(c);
            double diff = Math.abs(ic_a - ic_b);

            r += "\n" + c + "\t" + label + "\tic_a:" + ic_a + "\tic_b:" + ic_b + "\tdiff:" + diff;

        }
        return r;
    }

    private static void computeDiffDistrib(List<URI> listURIs, double[] ic1, double[] ic2, int nbClasses) {

        Map<Integer, Set<URI>> statClasses = new HashMap(nbClasses);

        for (int i = 0; i < ic1.length; i++) {

            int diff = (int) (Math.abs(ic1[i] - ic2[i]) * nbClasses);

            int classId = diff == nbClasses ? nbClasses - 1 : diff;

            if (statClasses.get(classId) == null) {
                statClasses.put(classId, new HashSet<URI>());
            }
            statClasses.get(classId).add(listURIs.get(i));
        }

        for (int i = 0; i < nbClasses; i++) {
            double sizeClass = statClasses.containsKey(i) ? statClasses.get(i).size() : 0;
            System.out.println("Class [" + ((double) i / nbClasses) + ";" + ((double) (i + 1) / nbClasses) + "]\t" + i + "\t" + sizeClass + "\t" + sizeClass / (double) ic1.length * 100.0);
            
        }

    }

    private static int[] getKMostVariant(int k, double[] tabk, double[] tabl) {

        Map<Integer, Double> diff = new HashMap(tabk.length);

        for (int i = 0; i < tabk.length; i++) {

            diff.put(i, Math.abs(tabk[i] - tabl[i]));

        }

        int[] k_mostVariants = new int[k];
        int i = 0;
        Map<Integer, Double> sortedDiff = MapUtils.sortByValueDecreasing(diff);

        for (Map.Entry<Integer, Double> e : sortedDiff.entrySet()) {
            k_mostVariants[i] = e.getKey();
            i++;
            if (i == k) {
                break;
            }
        }

        return k_mostVariants;
    }

}
