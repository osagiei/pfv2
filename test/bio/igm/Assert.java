package bio.igm;

import java.util.ArrayList;
import java.util.List;

/**
 * Minimal assertion helper, so the test suite needs no external dependency.
 */
public final class Assert {

    private static final List<String> FAILURES = new ArrayList<>();
    private static int checks;

    private Assert() {
    }

    public static void isTrue(String what, boolean condition) {
        checks++;
        if (!condition) {
            FAILURES.add(what + ": expected true");
        }
    }

    public static void isFalse(String what, boolean condition) {
        checks++;
        if (condition) {
            FAILURES.add(what + ": expected false");
        }
    }

    public static void equals(String what, Object expected, Object actual) {
        checks++;
        if (expected == null ? actual != null : !expected.equals(actual)) {
            FAILURES.add(what + ": expected <" + expected + "> but was <" + actual + ">");
        }
    }

    public static void throwsException(String what, Class<? extends Throwable> type, Runnable body) {
        checks++;
        try {
            body.run();
            FAILURES.add(what + ": expected " + type.getSimpleName() + " but nothing was thrown");
        } catch (Throwable t) {
            if (!type.isInstance(t)) {
                FAILURES.add(what + ": expected " + type.getSimpleName() + " but got " + t);
            }
        }
    }

    public static int checks() {
        return checks;
    }

    public static List<String> failures() {
        return FAILURES;
    }
}
