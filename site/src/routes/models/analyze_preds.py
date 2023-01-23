from datetime import datetime

from matbench_discovery.data import PRED_FILENAMES, load_df_wbm_with_preds

today = f"{datetime.now():%Y-%m-%d}"
models = list(PRED_FILENAMES)

df = load_df_wbm_with_preds(models)

df_out = df[models].isna().sum().to_frame(name="missing_preds")
df_out["test_set_size"] = len(df)


df_out.to_json(f"{today}-pred-analysis.json", orient="index")
