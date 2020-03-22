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
public class convertNetData_To_MatrixData_HeterNet {
    public static void main(String[] args) throws Exception{

        
        //Drug-Disease
        //Input
        String BiNet_FileName = "/Users/admin/Data/UsableData/BipartiteNets/DrugDiseaseNet_PREDICT.txt";

        //Output
        String AssocMat12_FileName = "/Users/admin/Manuscripts/52 HDR - Cytoscape App/Data/UnMatched/DrugDiseaseMat_PREDICT.txt";
        String AssocMat21_FileName = "/Users/admin/Manuscripts/52 HDR - Cytoscape App/Data/UnMatched/DiseaseDrugMat_PREDICT.txt";
        
        
        
        int i,j;
        
        //Read BipartiteNet
        BufferedReader brBiNet = new BufferedReader(new FileReader(BiNet_FileName));
        Map<String, Map<String, Double>> ID12ID2WeightMap = new TreeMap<>();
        Map<String, Map<String, Double>> ID22ID1WeightMap = new TreeMap<>();
        
        Set<String> AssocID2Set =  new TreeSet<>();
        Set<String> AssocID1Set =  new TreeSet<>();
        
        String strBiNet = "";
        
        while ((strBiNet = brBiNet.readLine()) != null) {
            String[] s = strBiNet.split("\t");
            String ID1 = s[0].trim();
            String ID2 = s[2].trim();
            
            AssocID1Set.add(ID1);
            AssocID2Set.add(ID2);
            
            if(ID12ID2WeightMap.containsKey(ID1)){
                ID12ID2WeightMap.get(ID1).put(ID2, 1.0);
            }else{
                Map<String, Double> ID2WeightMap = new TreeMap<>();
                ID2WeightMap.put(ID2, 1.0);
                ID12ID2WeightMap.put(ID1, ID2WeightMap);
            }
                
            if(ID22ID1WeightMap.containsKey(ID2)){
                ID22ID1WeightMap.get(ID2).put(ID1, 1.0);
            }else{
                Map<String, Double> ID1WeightMap = new TreeMap<>();
                ID1WeightMap.put(ID1, 1.0);
                ID22ID1WeightMap.put(ID2, ID1WeightMap);
            }
                        
        }
        brBiNet.close();
        System.out.println("ID12ID2WeightMap.size(): " + ID12ID2WeightMap.size());
        System.out.println("ID22ID1WeightMap.size(): " + ID22ID1WeightMap.size());
        
        System.out.println("AssocID1Set.size(): " + AssocID1Set.size() + "\t" + AssocID1Set.toString());
        System.out.println("AssocID2Set.size(): " + AssocID2Set.size() + "\t" + AssocID2Set.toString());
        
        List<String> AssocID1List = new ArrayList<>();
        AssocID1List.addAll(AssocID1Set);
        
        List<String> AssocID2List = new ArrayList<>();
        AssocID2List.addAll(AssocID2Set);
        
        
        
        
        //Construct AssocMat12_FileName
        PrintWriter pwAssocMat12 = new PrintWriter(new FileOutputStream(AssocMat12_FileName), true);
        //pwDrDiAssoc.print("Dr/Di");
        for(i=0;i<AssocID2List.size();i++){
            String ID2 = AssocID2List.get(i);
            pwAssocMat12.print("\t" + ID2);
        }
        pwAssocMat12.println();
        for(i=0;i<AssocID1List.size();i++){
            String ID1 = AssocID1List.get(i);
            pwAssocMat12.print(ID1);
            for(j=0;j<AssocID2List.size();j++){
                String ID2 = AssocID2List.get(j);
                int ina =0;
                if(ID12ID2WeightMap.get(ID1).containsKey(ID2)){
                    ina =1;
                }
                pwAssocMat12.print("\t" + ina);
                
            }
            pwAssocMat12.println();
        }
        pwAssocMat12.close();
        
        
        //Construct AssocMat21_FileName
        PrintWriter pwAssocMat21 = new PrintWriter(new FileOutputStream(AssocMat21_FileName), true);
        //pwDrDiAssoc.print("Di/Dr");
        for(i=0;i<AssocID1List.size();i++){
            String ID1 = AssocID1List.get(i);
            pwAssocMat21.print("\t" + ID1);
        }
        pwAssocMat21.println();
        for(i=0;i<AssocID2List.size();i++){
            String ID2 = AssocID2List.get(i);
            pwAssocMat21.print(ID2);
            for(j=0;j<AssocID1List.size();j++){
                String ID1 = AssocID1List.get(j);
                int ina =0;
                if(ID22ID1WeightMap.get(ID2).containsKey(ID1)){
                    ina =1;
                }
                pwAssocMat21.print("\t" + ina);
                
            }
            pwAssocMat21.println();
        }
        pwAssocMat21.close();
    }
}
