from flask_wtf import FlaskForm
from wtforms import SubmitField, TextField, FormField, FileField, IntegerField, FieldList, SelectField, FloatField, \
    BooleanField
from wtforms.validators import InputRequired, ValidationError, StopValidation, AnyOf, Regexp, NumberRange, DataRequired
from wtforms import StringField


class CopulaFactorForm(FlaskForm):
    gibbs_sampling_n = IntegerField("Gibbs sampling", validators=[InputRequired()], default=1000, description="")
    gibbs_burn_in_n = IntegerField("Gibbs burn in", validators=[InputRequired()], default=500, description="")
    gibbs_first_random_seed_n = IntegerField("Gibbs first random seed", validators=[InputRequired()], default=365,
                                             description="")
    gibbs_random_seed_update_parameter_n = IntegerField("Gibbs random seed update parameter",
                                                        validators=[InputRequired()], default=10, description="")


class EdgeWeightForm(FlaskForm):
    bootstrap_n = IntegerField("Bootstrap", validators=[InputRequired()], default=10, description="")
    bootstrap_first_random_seed_n = IntegerField("Bootstrap first random seed", validators=[InputRequired()],
                                                 default=100,
                                                 description="")
    bootstrap_random_seed_update_parameter_n = IntegerField("Bootstrap random seed update parameter",
                                                            validators=[InputRequired()], default=10, description="")


class PcAlgorithmForm(FlaskForm):
    causal_discovery_observation_n = IntegerField("Causal discovery observation", validators=[InputRequired()],
                                                  default=0, description="")
    indepTest = SelectField("indepTest",
                            choices=[('dsepTest', 'dsepTest'), ('disCItest', 'disCItest'), ('binCItest', 'binCItest'),
                                     ('gaussCItest', 'gaussCItest')],

                            default='gaussCItest')

    alpha = FloatField("Alpha", validators=[InputRequired()], default=0.05, description="")
    numCores = IntegerField("Number cores", validators=[InputRequired()], default=1, description="")
    verbose = SelectField("Verbose", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                           coerce=lambda x: x == 'True',
                          default='False')
    # fixedGaps = StringField("FixedGaps", validators=[InputRequired()], default='NULL')
    # fixedEdges = StringField("FixedEdges", validators=[InputRequired()], default='NULL')
    NAdelete = SelectField("NAdelete", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                           coerce=lambda x: x == 'True', default='True')
    m_max = StringField("M.max",  validators=[InputRequired(), Regexp(r'^(\d+|[iI][Nn][Ff])$')], default='Inf')

    u2pd = SelectField("u2pd", choices=[('relaxed', 'relaxed'), ('rand', 'rand'), ('retry', 'retry')], default='relaxed')
    skel_method = SelectField("skel_method",
                              choices=[('stable', 'stable'), ('original', 'original'), ('stable.fast', 'stable.fast')],
                              default='stable')

    conservative = SelectField("Conservative", choices=[(True, 'True'), (False, 'False')],validators=[InputRequired()],
                           coerce=lambda x: x == 'True',  default='False')
    maj_rule = SelectField("Maj rule", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                           coerce=lambda x: x == 'True', default='True')
    solve_confl = SelectField("Solve confl", choices=[(True, 'True'), (False, 'False')], validators=[InputRequired()],
                              coerce=lambda x: x == 'True', default='True')


class PlotAndDisplayForm(FlaskForm):
    core_plot_title_str = StringField("Core plot title", validators=[InputRequired()], default='Fake model')


class GeneralParameterForm(FlaskForm):
    copula_factor = FormField(CopulaFactorForm)
    edge_weight = FormField(EdgeWeightForm)
    pc_algorithm = FormField(PcAlgorithmForm)
    plot_and_display = FormField(PlotAndDisplayForm)
