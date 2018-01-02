package main.hmm;

class Result {
    private final double probability;
    private final int argument;

    public Result(double probability, int argument) {
        this.probability = probability;
        this.argument = argument;
    }

    public double getProbability() {
        return probability;
    }

    public int getArgument() {
        return argument;
    }
}