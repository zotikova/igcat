# 0.
from sklearn import metrics
import pandas as pd

def merge(l1, l2):
    combined = []
    c1 = 0
    c2 = 0
    n1 = l1[c1]
    n2 = l2[c2]
    len1 = len(l1)
    len2 = len(l2)
    while (c1<len1 and c2<len2):
        if n1 < n2:
            combined.append(n1)
            c1+=1
            if (c1 < len1):
                n1 = l1[c1]
        else:
            combined.append(n2)
            c2+=1
            if (c2 < len2):
                n2 = l2[c2]
    combined.extend(l1[c1:])
    combined.extend(l2[c2:])
    return combined

def mergesort(a):
    n = len(a)
    if n==0:
        return a
    if n ==1:
        return a
    if n > 1:
        return merge(mergesort(a[:int(n/2)]), mergesort(a[int(n/2):]))

def get_values(p_fr1, n_fr1, fr11, tp, fp, fn, tn, i):
    if p_fr1 == list(range(fr11[i - 1], fr11[i] + 1)):
        tp.append(len(p_fr1))
        tn.append(len(n_fr1))
        fn.append(0)
        fp.append(0)
    elif p_fr1 < list(range(fr11[i - 1], fr11[i] + 1)):
        tp.append(len(set(p_fr1) & set(range(fr11[i - 1], fr11[i] + 1))))
        tn.append(len(set(n_fr1) & (set(range(1, 2002)) ^ set(list(range(fr11[i - 1], fr11[i] + 1))))))
        fn.append(len(list(set(list(range(fr11[i - 1], fr11[i] + 1))) ^ set(p_fr1))))
        fp.append(0)
    elif p_fr1 > list(range(fr11[i - 1], fr11[i] + 1)):
        tp.append(len(set(p_fr1) & set(range(fr11[i - 1], fr11[i] + 1))))
        fp.append(len(set(p_fr1) ^ set(range(fr11[i - 1], fr11[i] + 1))))
        tn.append(len(set(n_fr1) & (set(range(1, 2002)) ^ set(list(range(fr11[i - 1], fr11[i] + 1))))))
        fn.append(0)

# values for begining aa
def beg_metr(assess, initial):
    tp, fp, fn = 0, 0, 0
    for i in range(len(assess)):
        if i%2 == 0 and assess[i] == initial[i]:
            tp = tp + 1
        elif i%2 == 0 and assess[i] < initial[i]:
            fp = fp + 1
        elif i%2 == 0 and assess[i] > initial[i]:
            fn = fn + 1
    return tp, fp, fn

# values for ending aa
def end_metr(assess, initial):
    tp, fp, fn = 0, 0, 0
    for i in range(len(assess)):
        if i%2 != 0 and assess[i] == initial[i]:
            tp = tp + 1
        elif i%2 != 0 and assess[i] > initial[i]:
            fp = fp + 1
        elif i%2 != 0 and assess[i] < initial[i]:
            fn = fn + 1
    return tp, fp, fn

# metrices
def precision(tp, fp):
    return tp/(tp+fp)
def recall(tp, fn):
    return tp/(tp+fn)
def specificity(TN, FP):
    return TN/(TN+FP)
def F1_score(Precision, Recall):
    return 2*(Precision*Recall)/(Precision+Recall)
def accuracy(TP, TN, FP, FN):
    return (TP+TN)/(TP+TN+FP+FN)
def FDR(fp, tp):
    return fp/(fp+tp)

df1 = pd.DataFrame({
        'Fragment': [],
        'TP': [],
        'FP': [],
        'FN': [],
        'Precision': [],
        'Sensitivity': [],
        'F1-score': []})

# Fragment 1
def calc_metr(assess, initial, text1, text2, df1): # text1="Fr1", text2="beginning"
    if text2 == 'beginning':
        pr01 = '**** {} {} ****'.format(text1, text2)
        f1_beg = beg_metr(assess, initial)
    elif text2 == 'ending':
        pr02 = '**** {} {} ****'.format(text1, text2)
        f1_beg = end_metr(assess, initial)
    pr1 = 'TP = {}, FP = {}, FN = {}'.format(f1_beg[0], f1_beg[1], f1_beg[2])
    prec_f1_beg = precision(f1_beg[0], f1_beg[1])
    pr2 = 'Precision for {} {} is {}'.format(text1, text2, prec_f1_beg)
    recall_f1_beg = recall(f1_beg[0], f1_beg[2])
    pr3 = 'Sensitivity for {} {} is {}'.format(text1, text2, recall_f1_beg)
    F1_f1_beg = F1_score(prec_f1_beg, recall_f1_beg)
    pr4 = 'F1-score for {} {} is {} \n'.format(text1, text2, F1_f1_beg)
    df2 = pd.DataFrame({
        'Fragment': ['{} {}'.format(text1, text2)],
        'TP': [f1_beg[0]],
        'FP': [f1_beg[1]],
        'FN': [f1_beg[2]],
        'Precision': [prec_f1_beg],
        'Sensitivity': [recall_f1_beg],
        'F1-score': [F1_f1_beg]})
    df1 = df1.append(df2, ignore_index=True)
    #df1 = pd.concat(df1, df2)
    return df1

# 1. To reassign lens to extended seqs
def get_values2(file_to_compare, file_reference, fragment):
    with open(file_to_compare, 'r') as f1, open(file_reference, 'r') as f2:
        f1, f2 = f1.readlines(), f2.readlines()
        f1, f2 = mergesort(f1), mergesort(f2)
        fr11, fr21, fr31, cdr11, cdr21, cdr31 = [], [], [], [], [], []
        fr12, fr22, fr32, cdr12, cdr22, cdr32 = [], [], [], [], [], []
        for _ in f2:
            _ = _.split('\t')
            fr12 += (int(_[1]), int(_[2]))
            fr22 += (int(_[5]), int(_[6]))
            fr32 += (int(_[9]), int(_[10]))
            cdr12 += (int(_[3]), int(_[4]))
            cdr22 += (int(_[7]), int(_[8]))
            cdr32 += (int(_[11]), int(_[12]))
        for _ in f1:
            _ = _.split('\t')
            fr11 += (int(_[1]), int(_[2]))
            fr21 += (int(_[5]), int(_[6]))
            fr31 += (int(_[9]), int(_[10]))
            cdr11 += (int(_[3]), int(_[4]))
            cdr21 += (int(_[7]), int(_[8]))
            cdr31 += (int(_[11]), int(_[12]))

    tp, fp, fn, tn = [], [], [], []
    for i in range(len(fr12)):
        if i % 2 != 0:
            frag_len = list(range(1, cdr32[i] + 1))
            p_fr1 = list(range(fr12[i - 1], fr12[i] + 1))
            n_fr1 = set(frag_len) - set(p_fr1)

            p_cdr1 = list(range(cdr12[i - 1], cdr12[i] + 1))
            n_cdr1 = set(frag_len) - set(p_cdr1)

            p_fr2 = list(range(fr22[i - 1], fr22[i] + 1))
            n_fr2 = set(frag_len) - set(p_fr2)

            p_cdr2 = list(range(cdr12[i - 1], cdr12[i] + 1))
            n_cdr2 = set(frag_len) - set(p_cdr2)

            p_fr3 = list(range(fr32[i - 1], fr32[i] + 1))
            n_fr3 = set(frag_len) - set(p_fr3)

            p_cdr3 = list(range(cdr32[i - 1], cdr32[i] + 1))
            n_cdr3 = set(frag_len) - set(p_cdr3)
            # print(list(n_cdr3))

            if fragment == 'fr1':
                val = get_values(p_fr1, n_fr1, fr11, tp, fp, fn, tn, i)
            elif fragment == 'fr2':
                val = get_values(p_fr2, n_fr2, fr21, tp, fp, fn, tn, i)
            elif fragment == 'fr3':
                val = get_values(p_fr3, n_fr3, fr31, tp, fp, fn, tn, i)
            elif fragment == 'cdr1':
                val = get_values(p_cdr1, n_cdr1, cdr11, tp, fp, fn, tn, i)
            elif fragment == 'cdr2':
                val = get_values(p_cdr2, n_cdr2, cdr21, tp, fp, fn, tn, i)
            else:
                val = get_values(p_cdr3, n_cdr3, cdr31, tp, fp, fn, tn, i)
    return sum(tp), sum(tn), sum(fp), sum(fn)

# 2. To calculate metrics for particular set of statistical parameters and particular antibody region
def metrices(ros_values, k):
    prec_ros_seq = precision(ros_values[k][0], ros_values[k][2])
    rec_ros_seq = recall(ros_values[k][0], ros_values[k][3])
    spec_ros_sec = specificity(ros_values[k][1], ros_values[k][2])
    F1_ros_sec = F1_score(prec_ros_seq, rec_ros_seq)
    acc_ros_sec = accuracy(ros_values[k][0], ros_values[k][1], ros_values[k][2], ros_values[k][3])
    FDR_ros_sec = FDR(ros_values[k][2], ros_values[k][0])
    return prec_ros_seq, rec_ros_seq, spec_ros_sec, F1_ros_sec, acc_ros_sec, FDR_ros_sec

# 3. Printing chart
def get_chart(pl1, pl2, ind1, ind2, title, ind3):
    df_calc = pd.DataFrame({
        'group': [ind1, ind2, ind3],
        'Precision': [pl1[0], pl2[0], 1.0],
        'Sensitivity': [pl1[1], pl2[1], 1.0],
        'Specificity': [pl1[2], pl2[2], 1.0],
        'F1-score': [pl1[3], pl2[3], 1.0],
        'Accuracy': [pl1[4], pl2[4], 1.0],
        'FDR': [pl1[5], pl2[5], 0.0]
    })

    # ------- PART 1: Create background
    from matplotlib.pyplot import figure
    figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')

    # number of variable
    categories = list(df_calc)[1:]
    N = len(categories)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    ax = plt.subplot(111, polar=True)

    # If you want the first axis to be on top:
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, size=11)

    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([0.2, 0.4, 0.6, 0.8, 1], ["0.2", "0.4", "0.6", "0.8", "1"], color="grey", size=9)
    plt.ylim(0, 1)

    # ------- PART 2: Add plots

    # Plot each individual = each line of the data
    # I don't do a loop, because plotting more than 3 groups makes the chart unreadable

    # Ind1
    values = df_calc.loc[0].drop('group').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, 'yellow', linewidth=1, linestyle='solid', label=ind1)  # orchid
    ax.fill(angles, values, 'yellow', alpha=0.2)

    # Ind2
    values = df_calc.loc[1].drop('group').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, 'purple', linewidth=1, linestyle='solid', label=ind2)  # lime
    ax.fill(angles, values, 'purple', alpha=0.2)

    # Ind3
    values = df_calc.loc[2].drop('group').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, 'black', linewidth=2, linestyle='--', label=ind3)
    ax.fill(angles, values, 'black', alpha=0)

    # Add legend
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

    # Figure title
    plt.title(title, fontsize=14, fontweight='bold')

    df_calc.to_csv(r'/home/sola/df_Fr2.csv')

# 4. Application

ros_values_fr3 = get_values2('/home/sola/heavy_1000_rosss.marking', '/home/sola/vh.kabat', 'fr3')

frags_to_anal = ['fr1', 'fr2', 'fr3', 'cdr1', 'cdr2', 'cdr3']

ros_values = []
for i in frags_to_anal:
    ros_values.append(get_values2('/home/sola/heavy_1000_rosss.marking', '/home/sola/heavy_1000.marking', i))
print(ros_values)

ros_metr = [prec_ros_seq, rec_ros_seq, spec_ros_sec, F1_ros_sec, acc_ros_sec, FDR_ros_sec]

igcat_values = []
for i in frags_to_anal:
    igcat_values.append(get_values2('/home/sola/heavy_1000_ini.marking', '/home/sola/heavy_1000.marking', i))
print(igcat_values)

igcat_metr = metrices(igcat_values, 1)
print(igcat_metr)

get_chart(igcat_metr, ros_metr, 'Fr2 IgCAT', 'Fr2 RosettaAntibody', 'Subunit Fr2', 'Fr2 IgBLAST')

get_chart(igcat_metr, ros_metr, 'Fr3 IgCAT', 'Fr3 RosettaAntibody', 'Subunit Fr3', 'Fr3 IgBLAST')

mini_ros = [(18150, 59770, 1178, 0), (6947, 50244, 28, 4508), (18222, 38808, 118, 4123)]
glob_ros = map(sum, zip(*mini_ros))
glob_ros2 = []
for i in glob_ros:
    glob_ros2.append(i)
glob_ros2 = [glob_ros2]
glob_ros3 = metrices(glob_ros2, 0)

get_chart(glob_cat3, glob_ros3, 'IgCAT', 'RosettaAntibody', 'All subunits', 'IgBLAST')