open module meltt {
    requires beast.pkgmgmt;
    requires beast.base;
    requires org.apache.commons.statistics.distribution;
    requires java.xml;
    requires javafx.base;
    requires java.desktop;

    exports meltt;

    provides beast.base.core.BEASTInterface with
            meltt.LTTLikelihood,
            meltt.LTT,
            meltt.CoalescentLTTPrior,
            meltt.AdjacentDeltaExchangeOperator;
}
