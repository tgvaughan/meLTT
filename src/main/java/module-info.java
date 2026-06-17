open module toposmc {
    requires beast.pkgmgmt;
    requires beast.base;
    requires org.apache.commons.statistics.distribution;
    requires java.xml;

    exports toposmc;

     provides beast.base.core.BEASTInterface with
         toposmc.LTTLikelihood;
}
