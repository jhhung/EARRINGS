def cal_matrix(TP, TN, FP, FN):
    if TN == 0:
        print("Add 1 to TN as smoothing.")
        TN += 1
    if TP == 0:
        print("Add 1 to TP as smoothing.")
        TP += 1
    if FN == 0:
        print("Add 1 to FN as smoothing.")
        FN += 1
    if FP == 0:
        print("Add 1 to FP as smoothing.")
        FP += 1

    ACC = (TP + TN) / (TP + TN + FP + FN)

    TPR = TP / (TP + FN)
    TNR = TN / (TN + FP)
    PPV = TP / (TP + FP)
    NPV = TN / (TN + FN)

    FNR = 1 - TPR
    FPR = 1 - TNR
    FDR = 1 - PPV
    FOR = 1 - NPV

    MCC = pow(PPV*TPR*TNR*NPV, 0.5) - pow(FDR*FNR*FPR*FOR, 0.5)

    matrix = [TP, TN, FP, FN]

    print("[TP, TN, FP, FN]:", TP, "," ,TN, ",", FP, ",", FN)
    print("ACC:", ACC)
    print("SEN:", TPR)
    print("SPC:", TNR)
    print("PPV:", PPV)
    print("MCC:", MCC)

    return ACC, TPR, TNR, PPV, MCC
