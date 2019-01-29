def is_methylated(base):
    """Returns True if the base denotes methylated cytosine.

    :param int base: Binary representation of methylation state of a cytosine.

    :returns: True if base == 1 else False
    """
    return base == 1

class HMMModel:
    """DNMT1-like HMM model for in silico proofreading.
    """
    def __init__(self, template_pattern, e_m=0.010, e_d=0.001, e_b=0.01, p=0.960, q=0.96):
        self.e = self._emission_probabilities(template_pattern, e_m, e_d, e_b)  # Emission probabilities.
        self.e_m, self.e_d, self.e_b = e_m, e_d, e_b
        self.init = [0.5, 0.5]  # Initial probabilities.
        self.p = p  # Processivity.
        self.q = q  # Recruitment efficiency.

    def _emission_probabilities(self, template_pattern, e_m, e_d, e_b):
        emission = []
        for base in template_pattern:
            if is_methylated(base):
                # [[e_dm(u), e_dm(m)], [e_am(u), e_am(m)]]
                emission.append([[1-e_d, e_d], [e_m, 1-e_m]])
            else:
                # [[e_du(u), e_du(m)], [e_au(u), e_au(m)]]
                emission.append([[1-e_b, e_b], [1-e_m, e_m]])

        return emission

    def proba(self, seq):
        """Returns the probability of the observed methylation pattern sequence given the template pattern.

        :param list seq: Binary representation of methylation pattern sequence.

        :returns: Probability of the observed methylation pattern based on the model.
        """
        prev = self.init
        p, q = self.p, self.q
        for i in range(len(seq)):
            curr = [self.e[i][0][seq[i]] * (prev[0] * (1 - q) + prev[1] * (1 - p)), self.e[i][1][seq[i]] * (prev[0] * q + prev[1] * p)]
            prev = curr

        s = sum(curr)

        prev = self.init
        for i in range(len(seq)-1, -1, -1):
            curr = [self.e[i][0][seq[i]] * (prev[0] * (1 - q) + prev[1] * (1 - p)), self.e[i][1][seq[i]] * (prev[0] * q + prev[1] * p)]
            prev = curr

        return (s + sum(curr)) / 2