package br.ufscar.hclga.classes;

/**
 *
 * @author Gean Trindade <gean.pereira@ufscar.br / geantrinpereira@gmail.com>
 */
public class Chronometer {

public long msBegin, msEnd, msTime;

public void start() {
    msTime = 0;
    msBegin = System.currentTimeMillis();
}

public void stop() {
    msEnd = System.currentTimeMillis();
    msTime += msEnd - msBegin;
}

public void resume() {
    msBegin = System.currentTimeMillis();
}

public long time() {
    return msTime;
}

public double stime() {
    return msTime / 1000.;
}

public double mtime() {
    return msTime / 60000.;
}

public double htime() {
    return msTime / 3600000.;
}

}
