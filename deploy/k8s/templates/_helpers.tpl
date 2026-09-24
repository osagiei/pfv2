{{- define "pfv2.name" -}}
{{- default .Chart.Name .Values.nameOverride | trunc 63 | trimSuffix "-" -}}
{{- end -}}

{{- define "pfv2.fullname" -}}
{{- printf "%s-%s" .Release.Name (include "pfv2.name" .) | trunc 52 | trimSuffix "-" -}}
{{- end -}}

{{- define "pfv2.labels" -}}
app.kubernetes.io/name: {{ include "pfv2.name" . }}
app.kubernetes.io/instance: {{ .Release.Name }}
app.kubernetes.io/version: {{ .Chart.AppVersion | quote }}
app.kubernetes.io/managed-by: {{ .Release.Service }}
{{- end -}}

{{- define "pfv2.serviceAccountName" -}}
{{- if .Values.serviceAccount.create -}}
{{- default (include "pfv2.fullname" .) .Values.serviceAccount.name -}}
{{- else -}}
{{- default "default" .Values.serviceAccount.name -}}
{{- end -}}
{{- end -}}

{{/* A digest pins what actually runs; fall back to the tag, then to the chart version. */}}
{{- define "pfv2.image" -}}
{{- if .Values.image.digest -}}
{{- printf "%s@%s" .Values.image.repository .Values.image.digest -}}
{{- else -}}
{{- printf "%s:%s" .Values.image.repository (default .Chart.AppVersion .Values.image.tag) -}}
{{- end -}}
{{- end -}}

{{/*
Job name for one sample. Sample ids commonly share a long prefix (SAMEA120815325 and
SAMEA120815326), so a truncated id is not unique; the hash suffix makes the name safe for
any id and any length, while the readable prefix keeps `kubectl get jobs` legible.
*/}}
{{- define "pfv2.jobName" -}}
{{- $prefix := include "pfv2.fullname" .root | trunc 24 | trimSuffix "-" -}}
{{- $slug := regexReplaceAll "[^a-z0-9-]" (.sample.id | lower) "-" | trunc 24 | trimSuffix "-" -}}
{{- printf "%s-%s-%s" $prefix $slug (.sample.id | sha256sum | trunc 8) -}}
{{- end -}}
