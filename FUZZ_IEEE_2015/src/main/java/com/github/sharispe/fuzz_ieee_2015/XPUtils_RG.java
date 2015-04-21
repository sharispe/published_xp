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
import slib.sml.sm.core.utils.SMconf;

/**
 *
 * @author Sébastien Harispe <sebastien.harispe@gmail.com>
 */
public class XPUtils_RG {

    public static void main(String[] args) throws Exception {

        boolean useOnlyNouns = true;
        boolean buildConceptRepresentation = false; // if the file Conf.file_concept_representation does not exists concept representations must be computed

        Set<BenchmarkEntry> benchmarkEntries = DataLoaderUtils.loadBenchmarkEntries_MC(Conf.benchmarkRG);

        URIFactory factory = URIFactoryMemory.getSingleton();
        G wordNet = DataLoaderUtils.loadWordNet(factory, Conf.dataloc, useOnlyNouns);
        DataLoaderUtils.loadAnnotations(factory, wordNet, Conf.annotations);

        Map<URI, String> labels = DataLoaderUtils.loadURILabels(factory, wordNet, Conf.dataloc, useOnlyNouns);

        FuzzIEEE_XPBrain computer = new FuzzIEEE_XPBrain(factory, wordNet, labels, Conf.file_concept_representation, buildConceptRepresentation, 1);

        SM_Engine engine = new SM_Engine(wordNet);

        ICconf icConf_sanchez = new IC_Conf_Topo(SMConstants.FLAG_ICI_SANCHEZ_2011);
        ICconf icConf_seco = new IC_Conf_Topo(SMConstants.FLAG_ICI_SECO_2004);
        ICconf icConf_resnik_i = new IC_Conf_Topo(SMConstants.FLAG_ICI_RESNIK_1995);
        ICconf icConf_resnik = new IC_Conf_Corpus(SMConstants.FLAG_IC_ANNOT_RESNIK_1995);

        SMconf conf_Lin_Sanchez = new SMconf("lin_sanchez", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998, icConf_sanchez);
        SMconf conf_Resnik_Sanchez = new SMconf("resnik_sanchez", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995, icConf_sanchez);

        SMconf conf_Lin_Seco = new SMconf("lin_seco", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998, icConf_seco);
        SMconf conf_Resnik_Seco = new SMconf("resnik_seco", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995, icConf_seco);

        SMconf conf_Lin_Resnik_i = new SMconf("lin_resnik_i", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998, icConf_resnik_i);
        SMconf conf_Resnik_Resnik_i = new SMconf("resnik_resnik_i", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995, icConf_resnik_i);

        SMconf conf_Lin_Resnik = new SMconf("lin_resnik", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998, icConf_resnik);
        SMconf conf_Resnik_Resnik = new SMconf("resnik_resnik", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995, icConf_resnik);

        List<SMconf> confs = new ArrayList();
        confs.add(conf_Lin_Sanchez);
        confs.add(conf_Resnik_Sanchez);
        confs.add(conf_Lin_Seco);
        confs.add(conf_Resnik_Seco);
        confs.add(conf_Lin_Resnik_i);
        confs.add(conf_Resnik_Resnik_i);
        confs.add(conf_Lin_Resnik);
        confs.add(conf_Resnik_Resnik);

        double[][] scoreMatrix = new double[confs.size()][benchmarkEntries.size()];
        double[] expectedScores = new double[benchmarkEntries.size()];

        double[] resnik_pignistic = new double[benchmarkEntries.size()];
        double[] lin_pignistic = new double[benchmarkEntries.size()];

        double[] resnik_plausibility = new double[benchmarkEntries.size()];
        double[] lin_plausibility = new double[benchmarkEntries.size()];

        double[] resnik_belief = new double[benchmarkEntries.size()];
        double[] lin_belief = new double[benchmarkEntries.size()];

        // initialize to negative value
        for (int i = 0; i < benchmarkEntries.size(); i++) {
            resnik_pignistic[i] = -1;
            lin_pignistic[i] = -1;

            resnik_plausibility[i] = -1;
            lin_plausibility[i] = -1;

            resnik_belief[i] = -1;
            lin_belief[i] = -1;
        }
        for (int i = 0; i < scoreMatrix.length; i++) {
            for (int j = 0; j < scoreMatrix[i].length; j++) {

                scoreMatrix[i][j] = -1;
            }
        }

        int i = 0;
        int j = 0;

        for (BenchmarkEntry b : benchmarkEntries) {

            URI uA = factory.getURI(Conf.prefixWordNet + b.label_A);
            URI uB = factory.getURI(Conf.prefixWordNet + b.label_B);

            i = 0;

            for (SMconf conf : confs) {

                double sim = engine.compare(conf, uA, uB);

                scoreMatrix[i][j] = sim;

                i++;
                System.out.println("\t" + conf.getLabel() + "\t" + b.label_A + "\t" + b.label_B + "\t"
                        + uA.getLocalName() + "\t" + uB.getLocalName() + "\texp:" + b.score
                        + "\tmeasure:" + sim);

            }

            // pignistic -----------------------------------------------------------------------------------------------
            double vlin_pignistic = XPUtils_MC_SIMLEX.computeLinMetric(uA, uB, engine, computer.getPignisticProbabilities());
            double vresnik_pignistic = XPUtils_MC_SIMLEX.computeResnikMetric(uA, uB, engine, computer.getPignisticProbabilities());

            // plausibility --------------------------------------------------------------------------------------------
            double vlin_plausibility = XPUtils_MC_SIMLEX.computeLinMetric(uA, uB, engine, computer.getPlausibility());
            double vresnik_plausibility = XPUtils_MC_SIMLEX.computeResnikMetric(uA, uB, engine, computer.getPlausibility());

            // belief -------------------------------------------------------------------------------------------------
            double vlin_belief = XPUtils_MC_SIMLEX.computeLinMetric(uA, uB, engine, computer.getBelief());
            double vresnik_belief = XPUtils_MC_SIMLEX.computeResnikMetric(uA, uB, engine, computer.getBelief());

            resnik_pignistic[j] = vresnik_pignistic;
            lin_pignistic[j] = vlin_pignistic;

            resnik_plausibility[j] = vresnik_plausibility;
            lin_plausibility[j] = vlin_plausibility;

            resnik_belief[j] = vresnik_belief;
            lin_belief[j] = vlin_belief;

            expectedScores[j] = b.score;
            j++;
        }

        CorrelationComputer corr = new CorrelationComputer();

        System.out.println(
                "----------------------------------------------------");
        System.out.println(
                "Pearson --------------------------------------------");
        System.out.println(
                "----------------------------------------------------");

        double corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, expectedScores);

        System.out.println(
                "checking coherence: " + corrPearson);

        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, resnik_pignistic);

        System.out.println(
                "\tresnik_pignistic\t" + corrPearson);
        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, lin_pignistic);

        System.out.println(
                "\tlin_pignistic\t" + corrPearson);

        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, resnik_belief);

        System.out.println(
                "\tresnik_belief\t" + corrPearson);
        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, lin_belief);

        System.out.println(
                "\tlin_belief\t" + corrPearson);

        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, resnik_plausibility);

        System.out.println(
                "\tresnik_plausibility\t" + corrPearson);
        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, lin_plausibility);

        System.out.println(
                "\tlin_plausibility\t" + corrPearson);
        i = 0;
        for (SMconf conf : confs) {

            double[] sims = Utils.extractColumn(i, scoreMatrix);
            corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.PEARSON, expectedScores, sims);

            System.out.println("\t" + conf.getLabel() + "\t" + corrPearson);
            i++;
        }

        System.out.println(
                "----------------------------------------------------");
        System.out.println(
                "Spearman --------------------------------------------");
        System.out.println(
                "----------------------------------------------------");

        double corrSpearman = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, expectedScores);

        System.out.println(
                "checking coherence: " + corrSpearman);

        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, resnik_pignistic);

        System.out.println(
                "\tresnik_pignistic\t" + corrPearson);
        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, lin_pignistic);

        System.out.println(
                "\tlin_pignistic\t" + corrPearson);

        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, resnik_belief);

        System.out.println(
                "\tresnik_belief\t" + corrPearson);
        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, lin_belief);

        System.out.println(
                "\tlin_belief\t" + corrPearson);

        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, resnik_plausibility);

        System.out.println(
                "\tresnik_plausibility\t" + corrPearson);
        corrPearson = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, lin_plausibility);

        System.out.println(
                "\tlin_plausibility\t" + corrPearson);

        i = 0;
        for (SMconf conf : confs) {

            double[] sims = Utils.extractColumn(i, scoreMatrix);
            corrSpearman = corr.computeCorrelation(CorrelationComputer.COR_METH.SPEARMAN, expectedScores, sims);

            System.out.println("\t" + conf.getLabel() + "\t" + corrSpearman);
            i++;
        }
    }

}
