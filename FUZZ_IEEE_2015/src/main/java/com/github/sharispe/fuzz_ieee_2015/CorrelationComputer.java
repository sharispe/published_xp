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

import au.com.bytecode.opencsv.CSVReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import slib.utils.ex.SLIB_Ex_Critic;

/**
 *
 * @author Sébastien Harispe <sebastien.harispe@gmail.com>
 */
public class CorrelationComputer {
    
    
    public enum COR_METH{
        PEARSON,SPEARMAN;
    }

    static Logger logger = LoggerFactory.getLogger(CorrelationComputer.class);
    Map<String, List<String>> columns;
    PearsonsCorrelation corPearson = new PearsonsCorrelation();
    SpearmansCorrelation corSpearman = new SpearmansCorrelation();

    public void loadData(String filepath, char separator) throws SLIB_Ex_Critic {

        try {

            columns = new HashMap();
            List<String> columnsNames = new ArrayList();

            //read file
            CSVReader reader = new CSVReader(new FileReader(filepath), separator);


            logger.info("Loading input from " + filepath);
            String[] nextLine = null;
            boolean header = true;

            while ((nextLine = reader.readNext()) != null) {

                if (header) {
                    header = false;
                    for (String s : nextLine) {
                        if (columns.containsKey(s)) {
                            throw new SLIB_Ex_Critic("Duplicate key: " + s + " in header ..." + s);
                        }
                        columnsNames.add(s);
                        columns.put(s, new ArrayList<String>());
                    }
                } else {
                    for (int i = 0; i < nextLine.length; i++) {
                        columns.get(columnsNames.get(i)).add(nextLine[i]);
                    }
                }
            }

        } catch (Exception e) {
            throw new SLIB_Ex_Critic(e.getMessage());
        }
    }
    
    public double computeCorrelation(COR_METH corMethod, String headerFlag1, String headerFlag2) throws SLIB_Ex_Critic{
        
        double[] vals1 = getColumnsValuesAsDouble(headerFlag1);
        double[] vals2 = getColumnsValuesAsDouble(headerFlag2);
        
        if(corMethod == COR_METH.PEARSON){
            return corPearson.correlation(vals1, vals2);
        }
        else if(corMethod == COR_METH.SPEARMAN){
            return corSpearman.correlation(vals1, vals2);
        }
        else{
            throw new SLIB_Ex_Critic("Error unknown correlation method "+corMethod);
        }
    }
    
    public double computeCorrelation(COR_METH corMethod, double[] vals1, double[] vals2) throws SLIB_Ex_Critic{
        
        if(corMethod == COR_METH.PEARSON){
            return corPearson.correlation(vals1, vals2);
        }
        else if(corMethod == COR_METH.SPEARMAN){
            return corSpearman.correlation(vals1, vals2);
        }
        else{
            throw new SLIB_Ex_Critic("Error unknown correlation method "+corMethod);
        }
    }

    private double[] getColumnsValuesAsDouble(String headerFlag) throws SLIB_Ex_Critic {
        
        if(!columns.containsKey(headerFlag)){
            throw new SLIB_Ex_Critic(headerFlag+" as no associated values");
        }
        
        List<String> valsAsStrings = columns.get(headerFlag);
        double[] vals = new double[valsAsStrings.size()];
        
        for (int i = 0; i < valsAsStrings.size(); i++) {
            
            vals[i] = Double.parseDouble(valsAsStrings.get(i));
        }
        return vals;
    }
    
    

}
