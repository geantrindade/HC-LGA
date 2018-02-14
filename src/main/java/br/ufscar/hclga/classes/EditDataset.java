package br.ufscar.hclga.classes;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import weka.core.Instances;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class EditDataset {

public static void main(String[] args) throws Exception {

    String pathDataset = leConfig()[0];
    String nameDataset = leConfig()[1];
    String trainOrTest = "train";
    String sl = leConfig()[2];
    int folds = Integer.valueOf(leConfig()[2]);

    for (int i = 1; i <= folds; i++) {
        String fullPathFile = pathDataset + nameDataset + i + trainOrTest + ".arff";
        Instances data = new Instances(new BufferedReader(new FileReader(fullPathFile)));
        data.setClassIndex(data.numAttributes() - 1);

//        Attribute att = data.attribute("id");
//
//        if (att != null) {
//            data.deleteAttributeAt(att.index());
//        }
//
//        att = data.attribute("total bases");
//
//        if (att != null) {
//            data.deleteAttributeAt(att.index());
//        }

        for (int j = 0; j < data.classAttribute().numValues(); j++) {
            String value = data.classAttribute().value(j);
            value = value.replace("R/", "");
            value = value.replace("R", "");
//            value = value.replace(".", "/");
            data.renameAttributeValue(data.classIndex(), j, value);
        }

        ArrayList<Object> objects = Collections.list(data.classAttribute().enumerateValues());
        ArrayList<String> classElements = new ArrayList<>(objects.size());

        for (Object object : objects) {
            classElements.add(Objects.toString(object, null));
        }

        for (int j = 0; j < classElements.size(); j++) {
            String element = classElements.get(j);
            int indexBar = element.indexOf("/");

            if (indexBar != -1) {
                ArrayList<String> separatedClasses = new ArrayList<>();

                while (indexBar != -1) {
                    String aux = element.substring(0, indexBar);
                    separatedClasses.add(aux);

                    indexBar = element.indexOf("/", indexBar + 1);
                }

                for (int k = 0; k < separatedClasses.size(); k++) {
                    String test = separatedClasses.get(k);

                    if (!classElements.contains(test)) {
                        classElements.add(test);
                    }
                }
            }
        }

        Collections.sort(classElements);
        String wholeValue = "";

        for (int j = 0; j < classElements.size(); j++) {
            wholeValue += classElements.get(j) + ",";
        }

        int wv = wholeValue.length() - 1;
        wholeValue = wholeValue.substring(0, wv);
        data.renameAttribute(data.classAttribute(), wholeValue);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fullPathFile))) {
            writer.write(data.toString());
            writer.newLine();
            writer.flush();
            writer.close();
        }
        
        if(i == folds && trainOrTest.equalsIgnoreCase("train")){
            trainOrTest = "test";
            i = 0;
        }
    }
}

public static String[] leConfig() {
    String regExp[] = {"pathDatasetsaveInDirectory =",
        "nameDataset =",
        "folds ="};

    Pattern comment = Pattern.compile("#");
    String[] result = new String[3];

    for (int i = 0; i < regExp.length; i++) {
        try {
            FileReader reader = new FileReader("C:\\Users\\gean_\\Dropbox\\posGrad\\GACerriMaven\\src\\main\\java\\hmc_ga\\editaDataset.txt");
//            FileReader reader = new FileReader("editaDataset.txt");
//            FileReader reader = new FileReader("/home/geantrindade/Dropbox/posGrad/GACerriMaven/src/main/java/hmc_ga/validacao.txt");
            BufferedReader buffReader = new BufferedReader(reader);
            Pattern pattern = Pattern.compile(regExp[i]);
            String line = null;

            while ((line = buffReader.readLine()) != null) {
                Matcher m = pattern.matcher(line);
                Matcher m1 = comment.matcher(line);

                if (m.find() && !m1.find()) {
                    String[] vectorLine = line.split(" = ");
                    result[i] = vectorLine[1];
                    break;
                }
            }

            buffReader.close();
            reader.close();

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    return result;
}
}
