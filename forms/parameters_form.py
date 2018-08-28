from flask_wtf import FlaskForm
from wtforms import SubmitField, TextField, FormField, FileField, IntegerField, FieldList, SelectField, FloatField, \
    HiddenField
from wtforms.validators import InputRequired, ValidationError, StopValidation, AnyOf, Regexp, NumberRange, DataRequired
from wtforms import StringField


class CopulaFactorForm(FlaskForm):
    gibbs_sampling_n = IntegerField("Number of Gibbs samples", validators=[InputRequired()], default=1000, description="")
    gibbs_burn_in_n = IntegerField("Number of Gibbs burn-in samples", validators=[InputRequired()], default=500, description="")
    gibbs_first_random_seed_n = IntegerField("Gibbs sampling first random seed", validators=[InputRequired()], default=365,
                                             description="")
    gibbs_random_seed_update_parameter_n = IntegerField("Gibbs sampling random seed update parameter",
                                                        validators=[InputRequired()], default=10, description="")


class EdgeWeightForm(FlaskForm):
    bootstrap_n = IntegerField("Number of bootstrap samples", validators=[InputRequired()], default=10, description="")
    bootstrap_first_random_seed_n = IntegerField("Bootstrapping first random seed", validators=[InputRequired()],
                                                 default=100,
                                                 description="")
    bootstrap_random_seed_update_parameter_n = IntegerField("Bootstrapping random seed update parameter",
                                                            validators=[InputRequired()], default=10, description="")


class PcAlgorithmForm(FlaskForm):
    # causal_discovery_observation_n : when “automatic” is selected,
    #  the `INI` file will have the line `causal_discovery_observation_n=0`
    causal_discovery_observation_n = StringField("Causal discovery observation number",
                                                 validators=[InputRequired(),
                                                             Regexp(r'^(\d+|[aA][uU][tT][oO][mM][aA][tT][iI][cC])$')],
                                                 default='automatic', description="")

    indepTest = SelectField("indepTest",
                            choices=[('dsepTest', 'dsepTest'), ('disCItest', 'disCItest'), ('binCItest', 'binCItest'),
                                     ('gaussCItest', 'gaussCItest')],

                            default='gaussCItest', description="Function for testing conditional independence.")

    alpha = FloatField("Alpha", validators=[InputRequired(), NumberRange(min=0, max=1)], default=0.05,
                       description="Significance level (number in"
                                   "(0,1) for the individual "
                                   "conditional independence "
                                   "tests.")
    numCores = IntegerField("Number cores", validators=[InputRequired()], default=1, description="Specifies the number"
                                                                                                 " of cores to be used "
                                                                                                 "for parallel "
                                                                                                 "estimation of "
                                                                                                 "skeleton.")
    verbose = SelectField("Verbose", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                          coerce=lambda x: x == 'True',
                          default='False', description="If TRUE, detailed output is provided.")

    fixedGaps = HiddenField("FixedGaps", validators=[InputRequired()], default='NULL')
    fixedEdges = HiddenField("FixedEdges", validators=[InputRequired()], default='NULL')

    NAdelete = SelectField("NAdelete", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                           coerce=lambda x: x == 'True', default='True',
                           description="If indepTest returns NA and this option is TRUE, the corresponding edge is "
                                       "deleted. If this option is FALSE, the edge is not deleted.")
    m_max = StringField("M.max", validators=[InputRequired(), Regexp(r'^(\d+|[iI][Nn][Ff])$')], default='Inf',
                        description="Maximal size of the conditioning sets that are considered in the conditional "
                                    "independence tests. (Inf = infinite)")

    u2pd = SelectField("u2pd", choices=[('relaxed', 'relaxed'), ('rand', 'rand'), ('retry', 'retry')],
                       default='relaxed',
                       description="String specifying the method for dealing with conflicting information when trying "
                                   "to orient edges.")
    skel_method = SelectField("skel_method",
                              choices=[('stable', 'stable'), ('original', 'original'), ('stable.fast', 'stable.fast')],
                              default='stable',
                              description="Character string specifying method; the default, 'stable' provides an "
                                          "order-independent skeleton.")

    conservative = SelectField("Conservative", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                               coerce=lambda x: x == 'True', default='False',
                               description="Logical indicating if the conservative PC is used.")
    maj_rule = SelectField("Maj rule", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                           coerce=lambda x: x == 'True', default='True',
                           description="Logical indicating that the triples shall be checked for ambiguity using a "
                                       "majority rule idea, which is less strict than the conservative PC algorithm.")
    solve_confl = SelectField("Solve confl", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                              coerce=lambda x: x == 'True', default='True',
                              description="If TRUE, the orientation of the v-structures and the orientation rules work "
                                          "with lists for candidate sets and allow bi-directed edges to resolve "
                                          "conflicting edge orientations.")


class PlotAndDisplayForm(FlaskForm):
    core_plot_title_str = StringField("Core plot title", validators=[InputRequired()], default='Result')


class GeneralParameterForm(FlaskForm):
    plot_and_display = FormField(PlotAndDisplayForm)
    pc_algorithm = FormField(PcAlgorithmForm)
    copula_factor_algorithm = FormField(CopulaFactorForm)
    edge_weight_algorithm = FormField(EdgeWeightForm)
