package org.broadinstitute.hellbender.engine.filters;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContext.Type;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;

import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;
import java.util.*;

/**
 * Tests for the various VariantContext variant filterSuffix predicates
 */
public class VariantFilterUnitTest extends GATKBaseTest {

    final Allele SnpRef = Allele.create("A", true);
    final Allele Snp = Allele.create("T");
    final Allele MnpRef = Allele.create("AAAAAAAAA", true);
    final Allele Mnp = Allele.create("CCCCCCCCC");

    VariantContext snpVC;
    VariantContext mnpVC;

    public VariantFilterUnitTest() throws FileNotFoundException {
        initializeHG19Reference();
        snpVC = createArtificialVC(
                "id1",
                new SimpleInterval("1", 42, 42),
                Arrays.asList(SnpRef, Snp)
        );

        mnpVC = createArtificialVC(
                "id3",
                new SimpleInterval("2", 2, 10),
                Arrays.asList(MnpRef, Mnp)
        );
    }

    /**
     * Create an artificial VariantContext
     *
     */
    private static VariantContext createArtificialVC(
            String id,
            SimpleInterval loc,
            List<Allele> alleles)
    {
        VariantContextBuilder vb = new VariantContextBuilder();
        vb.id(id);
        if (alleles != null)
            vb.alleles(alleles);
        if (loc != null)
            vb.loc(loc.getContig(), loc.getStart(), loc.getEnd());
        return vb.make();
    }

    @DataProvider(name="includeIDsVCs")
    public Object[][] includeIDsTestVCs() {

        return new Object[][]{
                { snpVC, new String[]{"id1"}, true },
                { snpVC, new String[]{"id1", "id2"}, true },
                { snpVC, new String[]{"noid"}, false },
        };
    }

    @Test(dataProvider="includeIDsVCs")
    public void testIncludeIDsVariantFilter(VariantContext vc, String[] incIDs, boolean expected) {
        Set<String> idSet = new LinkedHashSet<>();
        idSet.addAll(Arrays.asList(incIDs));
        VariantIDsVariantFilter iivf = new VariantIDsVariantFilter(idSet);
        Assert.assertTrue(iivf.test(vc) == expected);
    }

    @DataProvider(name="typeVCs")
    public Object[][] typeTestVCs() {

        return new Object[][]{
                { snpVC, new Type[]{Type.SNP}, true },
                { snpVC, new Type[]{Type.SNP, Type.MNP}, true },
                { snpVC, new Type[]{Type.INDEL}, false },
                { snpVC, new Type[]{Type.INDEL, Type.MIXED}, false },
                { mnpVC, new Type[]{Type.SNP}, false },
                { mnpVC, new Type[]{Type.SNP, Type.MNP}, true },
                { mnpVC, new Type[]{Type.INDEL}, false },
                { mnpVC, new Type[]{Type.INDEL, Type.MIXED}, false },
        };
    }

    @Test(dataProvider="typeVCs")
    public void testVariantTypeVariantFilter(VariantContext vc, Type[] types, boolean expected) {
        Set<Type> typesSet = new LinkedHashSet<>();
        typesSet.addAll(Arrays.asList(types));
        VariantTypesVariantFilter vtvf = new VariantTypesVariantFilter(typesSet);
        Assert.assertTrue(vtvf.test(vc) == expected);
    }

    @DataProvider
    public Object[][] getNotSymbolicOrSVFilterTestVCs() {
        return new Object[][]{
                {snpVC, true},
                {new VariantContextBuilder(snpVC).alleles(   // SYMBOLIC
                        Arrays.asList(
                                Allele.create("A", true),
                                Allele.create("<CNV>", false))).make(), false},
                {new VariantContextBuilder(snpVC).alleles(   // SV = > 150 indel
                        Arrays.asList(
                                Allele.create("A", true),
                                Allele.create("A" + Utils.dupChar('G', 151), false))).make(), false},
        };
    }

    @Test(dataProvider = "getNotSymbolicOrSVFilterTestVCs")
    public void testNotSVOrSymbolicVariantFilter(VariantContext vc, boolean expected) {
        Assert.assertEquals(VariantFilterLibrary.NOT_SV_OR_SYMBOLIC.test(vc), expected);
    }

    @DataProvider
    public Object[][] getPassesFiltersFilterTestVCs() {
        return new Object[][]{
                {snpVC, true},
                {new VariantContextBuilder(snpVC).unfiltered().make(), true},
                {new VariantContextBuilder(snpVC).passFilters().make(), true},
                {new VariantContextBuilder(snpVC).filter("unbelievable").make(), false},
        };
    }

    @Test(dataProvider = "getPassesFiltersFilterTestVCs")
    public void testPassesFiltersVariantFilter(VariantContext vc, boolean expected) {
        Assert.assertEquals(VariantFilterLibrary.PASSES_FILTERS.test(vc), expected);
    }

}
