package main.logger;

import java.text.SimpleDateFormat;

public class Log {

    private static SimpleDateFormat time_formatter = new SimpleDateFormat("HH:mm:ss.SSS");

    private static boolean printDebug = false;
    private static boolean printError = true;
    private static boolean printInfo = true;
    private static boolean addTime = false;

    public static void setPrintDebug(boolean printDebug) {
        Log.printDebug = printDebug;
    }

    public synchronized static void e(String text) {
        if (printError)
            System.err.print(time() + text);
    }

    public synchronized static void eLine(String text) {
        if (printError)
            System.err.println(time() + text);
    }

    public synchronized static void iLine() {
        if (printInfo)
            System.out.println();
    }

    public synchronized static void iLine(String text) {
        if (printInfo)
            System.out.println(time() + text);
    }

    public synchronized static void i(String text) {
        if (printInfo)
            System.out.print(text);
    }

    public synchronized static void i(int text) {
        i(String.valueOf(text));
    }

    public static void i(double v) {
        i(String.valueOf(v));
    }

    public synchronized static void dLine() {
        if (printDebug)
            System.out.println();
    }

    public synchronized static void dLine(char i) {
        if (printDebug)
            System.out.println(time() + i);
    }

    public synchronized static void dLine(int i) {
        if (printDebug)
            System.out.println(time() + i);
    }

    public synchronized static void dLine(String text) {
        if (printDebug)
            System.out.println(time() + text);
    }

    public synchronized static void d(String text) {
        if (printDebug)
            System.out.print(text);
    }

    public synchronized static void d(int text) {
        d(String.valueOf(text));
    }

    public synchronized static void d(double v) {
        d(String.valueOf(v));
    }

    private static String time() {
        if (!addTime)
            return "";
        return time_formatter.format(System.currentTimeMillis()) + "  ";

    }
}
