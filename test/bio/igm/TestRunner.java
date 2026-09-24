package bio.igm;

import bio.igm.entities.PutativeStructureTest;
import bio.igm.entities.ReadsTest;
import bio.igm.utils.discovery.DiscoveryTest;
import bio.igm.utils.filter.FilterTest;

/**
 * Runs the PFv2 test suite. Exits non-zero when anything fails.
 */
public final class TestRunner {

    public static void main(String[] args) throws Exception {
        ReadsTest.run();
        PutativeStructureTest.run();
        DiscoveryTest.run();
        FilterTest.run();

        if (Assert.failures().isEmpty()) {
            System.out.println("PASS  " + Assert.checks() + " checks");
            return;
        }

        System.out.println("FAIL  " + Assert.failures().size() + " of " + Assert.checks() + " checks");
        for (String failure : Assert.failures()) {
            System.out.println("  - " + failure);
        }
        System.exit(1);
    }
}
