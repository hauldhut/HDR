package dti;


import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Duc-Hau
 */
public class convertNetData_To_MatrixData_HomoNet {
    public static void main(String[] args) throws Exception{
        //PPI
        //Input
        //String HomoNet1_FileName = "/Users/admin/Manuscripts/71 RWR_MH/Data/FLN_Extracted_Normalized_Percent_1.0.txt";
        //String SimMat1_FileName = "/Users/admin/Manuscripts/71 RWR_MH/Data/Mat/FLNGeneSimMat.txt";
        
        //String HomoNet1_FileName = "/Users/admin/Manuscripts/71 RWR_MH/Data/PPI.txt";
        //String SimMat1_FileName = "/Users/admin/Manuscripts/71 RWR_MH/Data/Mat/GeneSimMat.txt";
        
        //String HomoNet1_FileName = "/Users/admin/Data/UsableData/Drug_Data/HomoNets/DrugSimNet_CHEM.txt";
        //String SimMat1_FileName = "/Users/admin/Manuscripts/52 HDR - Cytoscape App/Data/UnMatched/DrugSimMat_CHEM.txt";
        
        //String HomoNet1_FileName = "/Users/admin/Data/UsableData/Disease_Data/DiseaseNet/OMIM/DiseaseSimNet_OMIM.txt";
        //String SimMat1_FileName = "/Users/admin/Manuscripts/52 HDR - Cytoscape App/Data/UnMatched/DiseaseSimMat_OMIM.txt";
        
        //Shared gene-based
        String HomoNet1_FileName = "/Users/admin/Data/UsableData/HomoNets/EnhancerNet/enh2target-1.0.2-Tab-AdjList-0.05.txt";        
        String SimMat1_FileName = "/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Data/UnMatched/EnhSimMat_SharedGene.txt";
        
        //Seq-based
        //String HomoNet1_FileName = "/Users/admin/Data/Enhancer/DiseaseEnhancer/enh2disease-1.0.2.txt.bed.fasta_simmat_alf_input.tsv_withHeader.txt_List.txt";
        //String SimMat1_FileName = "/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Data/UnMatched/EnhSimMat_Seq.txt";
        
        //DOSimNet
        //String HomoNet1_FileName = "/Users/admin/Data/UsableData/LncRNA_Data/HomoNets/DiseaseNet/DOBasedOMIMEntitySimilarityNet.txt";
        //String SimMat1_FileName = "/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Data/UnMatched/DOSimMat.txt";
        
        int i,j;
        
        //Read HomoNet1
        BufferedReader brHomoNet1 = new BufferedReader(new FileReader(HomoNet1_FileName));
        Map<String, Map<String, Double>> HomoNet1Map = new TreeMap<>();
        Set<String> NetID1Set = new TreeSet<>();
        String strHomoNet1 = "";
        while ((strHomoNet1 = brHomoNet1.readLine()) != null) {
            //System.out.println(strHomoNet1);
            String[] s = strHomoNet1.split("\t");
            String ID1 = s[0].trim();
            String ID2 = s[2].trim();
            double Weight = Double.parseDouble(s[1].trim());
            
                    
            NetID1Set.add(ID1);
            NetID1Set.add(ID2);
            
            if(HomoNet1Map.containsKey(ID1)){
                HomoNet1Map.get(ID1).put(ID2, Weight);
            }else{
                Map<String, Double> ID2WeightMap = new TreeMap<>();
                ID2WeightMap.put(ID2, Weight);
                HomoNet1Map.put(ID1, ID2WeightMap);
            }
        }
        brHomoNet1.close();
        System.out.println("HomoNet1Map.size(): " + HomoNet1Map.size());
        System.out.println("NetID1Set.size(): " + NetID1Set.size());
        List<String> NetID1List = new ArrayList<>();
        NetID1List.addAll(NetID1Set);
        
        //Construct ID1 Sim Mat 
        PrintWriter pwSimMat1 = new PrintWriter(new FileOutputStream(SimMat1_FileName), true);
        for(i=0;i<NetID1List.size();i++){
            String ID1 = NetID1List.get(i);
            pwSimMat1.print("\t" + ID1);
        }
        pwSimMat1.println();
        
        
        for(i=0;i<NetID1List.size();i++){
            String ID1i = NetID1List.get(i);
            //System.out.println(ID1i + "\t" + i);
            pwSimMat1.print(ID1i);
            for(j=0;j<NetID1List.size();j++){
                String ID1j = NetID1List.get(j);
                double sim = 0;
                if(ID1i.compareTo(ID1j)==0) sim = 1.0;
                if(HomoNet1Map.containsKey(ID1i)){
                    if(HomoNet1Map.get(ID1i).containsKey(ID1j)){
                        sim = HomoNet1Map.get(ID1i).get(ID1j);
                    }
                }
                if(HomoNet1Map.containsKey(ID1j)){
                    if(HomoNet1Map.get(ID1j).containsKey(ID1i)){
                        sim = HomoNet1Map.get(ID1j).get(ID1i);
                    }
                }
                
                pwSimMat1.print("\t" + sim);

            }
            pwSimMat1.println();
        }
        pwSimMat1.close();
        
    }
}
