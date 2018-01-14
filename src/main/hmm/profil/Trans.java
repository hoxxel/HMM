package main.hmm.profil;

public class Trans {
    private final char startState;
    private final char endState;

    public Trans(char startState, char endState) {
        this.startState = startState;
        this.endState = endState;
    }

    public char getStartState() {
        return startState;
    }

    public char getEndState() {
        return endState;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null || obj.getClass() != Trans.class)
            return false;
        Trans trans = (Trans) obj;
        return equals(trans.getStartState(), trans.getEndState());
    }

    public boolean equals(char startState, char endState) {
        if (this.startState != startState || this.endState != endState)
            return false;
        return true;
    }

    @Override
    public String toString() {
        return startState + "" + endState;
    }
}
