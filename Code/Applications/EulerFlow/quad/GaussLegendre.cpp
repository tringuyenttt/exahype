#include "EulerFlow/quad/GaussLegendre.h"

const double exahype::quad::gaussLegendreMaxNodes = 10;

#if EXAHYPE_ORDER==0
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {1.0000000000000000};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.5000000000000000};
#elif EXAHYPE_ORDER==1
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.5000000000000000, 0.5000000000000000};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.2113248654051871, 0.7886751345948129};
#elif EXAHYPE_ORDER==2
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.2777777777777778, 0.4444444444444444, 0.2777777777777778};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.1127016653792583, 0.5000000000000000, 0.8872983346207417};
#elif EXAHYPE_ORDER==3
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.1739274225687273, 0.3260725774312732, 0.3260725774312732, 0.1739274225687273};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262};
#elif EXAHYPE_ORDER==4
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.1184634425280953, 0.2393143352496831, 0.2844444444444444, 0.2393143352496831, 0.1184634425280953};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.0469100770306680, 0.2307653449471585, 0.5000000000000000, 0.7692346550528415, 0.9530899229693319};
#elif EXAHYPE_ORDER==5
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.0856622461895851, 0.1803807865240703, 0.2339569672863455, 0.2339569672863455, 0.1803807865240703, 0.0856622461895851};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.0337652428984240, 0.1693953067668678, 0.3806904069584016, 0.6193095930415985, 0.8306046932331322, 0.9662347571015760};
#elif EXAHYPE_ORDER==6
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.0647424830844337, 0.1398526957446367, 0.1909150252525594, 0.2089795918367347, 0.1909150252525594, 0.1398526957446367, 0.0647424830844337};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.0254460438286208, 0.1292344072003028, 0.2970774243113014, 0.5000000000000000, 0.7029225756886985, 0.8707655927996972, 0.9745539561713792};
#elif EXAHYPE_ORDER==7
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.0506142681451723, 0.1111905172266916, 0.1568533229389448, 0.1813418916891810, 0.1813418916891810, 0.1568533229389448, 0.1111905172266916, 0.0506142681451723};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.0198550717512319, 0.1016667612931866, 0.2372337950418355, 0.4082826787521751, 0.5917173212478250, 0.7627662049581645, 0.8983332387068134, 0.9801449282487682};
#elif EXAHYPE_ORDER==8
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.0406371941807921, 0.0903240803474286, 0.1303053482014692, 0.1561735385200014, 0.1651196775006299, 0.1561735385200014, 0.1303053482014692, 0.0903240803474286, 0.0406371941807921};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.0159198802461870, 0.0819844463366821, 0.1933142836497048, 0.3378732882980955, 0.5000000000000000, 0.6621267117019045, 0.8066857163502952, 0.9180155536633179, 0.9840801197538130};
#elif EXAHYPE_ORDER==9
    const double exahype::quad::gaussLegendreWeights[EXAHYPE_ORDER+1] =
        {0.0333356721543462, 0.0747256745752621, 0.1095431812579907, 0.1346333596549986, 0.1477621123573765, 0.1477621123573765, 0.1346333596549986, 0.1095431812579907, 0.0747256745752621, 0.0333356721543462};
    const double exahype::quad::gaussLegendreNodes[EXAHYPE_ORDER+1]   =
        {0.0130467357414141, 0.0674683166555077, 0.1602952158504878, 0.2833023029353764, 0.4255628305091844, 0.5744371694908156, 0.7166976970646236, 0.8397047841495122, 0.9325316833444923, 0.9869532642585859};
#endif

